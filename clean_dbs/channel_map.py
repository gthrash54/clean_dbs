import io
import json
import re
import zipfile
from datetime import datetime
from pathlib import Path
import numpy as np

#PZ5 impedance checks never report above 536.87 kOhm (2138 hits at exactly this value across 408 local blocks, 2026-10-11), so it means "open"
PZ5_CEILING = 536.87


def channel_names(block, stream_name, n_channels):
    """
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    stream_name: the stream in result["streams"] (e.g. "dbs" for brainvision, "ecos" for tdt)
    n_channels: how many channels that stream has, used as a check
    Returns: list of names like "ecog_a1", "dbs1", "fdi_contra", or None if we can't say for sure
    Brainvision: the names the recorder wrote in the .vhdr (Ch1=ECOGA1,...), lowercased and with ecog written ecog_a1 / ecog_p1 everywhere.
      Checked against the hand-made <name>_channels.csv on all 38 recordings that have one (2026-10-11): 30 match, but
      ug3005 stage1 erpassay has A and P swapped in the csv, the stage2 erp files on ug3003/ug3006 leave out ch9 (FDI_ipsi),
      and the csv renames some emg (FDI -> fdi_contra), so the .vhdr is used (it's on all 151 recordings and matches what MNE reads).
    TDT: always None for now. The Synapse experiment saved in each block's .tin says ecos ch k = PZ5 input k (one 28 or 20 channel
      sub-amp), but nothing in it (or the downloaded experiments or the methods doc) says which strip or lead contact is plugged
      into which PZ5 input. That wiring has to come from the OR sheets / headbox adapter, so we don't guess.
    """
    block = Path(block)
    if block.suffix.lower() != ".vhdr" or stream_name != "dbs": #read_bv puts everything in one stream called dbs
        return None
    names = re.findall(r"^Ch\d+=([^,]*),", block.read_text(errors="ignore"), re.M)
    if len(names) != n_channels:
        return None
    #"ECOGA1", "ECOG_A1" and "ecog_a1" are all used for the same thing depending on the patient
    names = [re.sub(r"^ecog_?([ap])(\d+)$", r"ecog_\1\2", n.strip().lower()) for n in names]
    return names if len(set(names)) == len(names) else None


def connected_contacts(block):
    """
    Reads the impedance checks saved with a recording, so you can tell which contacts were actually plugged in.
    Input:
    block: a tdt block folder or a brainvision .vhdr file
    Returns: {stream name: dictionary of arrays}, one entry per channel per check:
      channel: channel number in that stream (0 = the reference)
      time: when the check ran, seconds from the start of the recording (brainvision checks happen before, so negative)
      frequency: test frequency in Hz (nan for brainvision, it doesn't say)
      impedance: kOhm (nan when the rig only wrote a word like "Out of Range!")
      connected: 1 = measured below the ceiling, 0 = at the ceiling / out of range / disconnected
    TDT: PZ5(1)p<n>.csv is PZ5 sub-amp n; the block's .tin says which store that sub-amp feeds (e.g. p1 -> ecos). If the channel counts
      don't line up (emg sub-amps take 8 inputs but send 4 channels) or there's no .tin, the key is pz5_p<n> instead of a store name.
      IZV10 csvs (stimulator outputs, not recording channels) come back as izv10_p<n>. -1 in the csv means not measured in that pass.
    Heads up: a few open channels read just under the ceiling (534 to 536 kOhm on uh3008), so look at the impedance, not just connected.
    """
    block = Path(block)
    if block.suffix.lower() == ".vhdr":
        return _bv_impedance(block)
    stores = _pz5_stores(block)
    out = {}
    for csv in sorted(block.glob("*.csv")):
        match = re.fullmatch(r"(PZ5|IZV10)\(\d+\)p(\d+)\.csv", csv.name)
        if match is None:
            continue
        lines = [line.split(",") for line in csv.read_text(errors="ignore").splitlines() if line.strip()]
        if len(lines) < 2: #header only = nobody ran a check during this block
            continue
        #columns R1..Rn are the channels, REF the reference; TIME/FREQUENCY/TARGET are per check
        columns = {name.split(" ")[0]:i for i, name in enumerate(lines[0])}
        chans = [(int(name[1:]), i) for name, i in columns.items() if re.fullmatch(r"R\d+", name)]
        if "REF" in columns:
            chans.append((0, columns["REF"]))
        rows = {"channel":[], "time":[], "frequency":[], "impedance":[], "connected":[]}
        for line in lines[1:]:
            for number, i in chans:
                value = float(line[i])
                if value < 0: #-1 = this channel wasn't measured in this pass
                    continue
                rows["channel"].append(number)
                rows["time"].append(float(line[columns["TIME"]]))
                rows["frequency"].append(float(line[columns["FREQUENCY"]]))
                rows["impedance"].append(value)
                rows["connected"].append(int(value < PZ5_CEILING))
        n_inputs = len([c for c in chans if c[0] > 0])
        key = f"{match.group(1).lower()}_p{match.group(2)}"
        if match.group(1) == "PZ5" and stores.get(int(match.group(2)), (None, 0))[1] == n_inputs:
            key = stores[int(match.group(2))][0]
        out[key] = _as_arrays(rows)
    return out


def _pz5_stores(block):
    """
    Opens the Synapse experiment saved inside the block's .tin (zip -> .synexpz zip -> .synexp json).
    Returns: {PZ5 sub-amp number: (store name, channel count)} for stores that take a whole sub-amp starting at channel 1, e.g. {1:("ecos", 28)}
    Found this way on all 408 local blocks (2026-10-11): ecogdbs/ecog/ecogamp is always sub-amp 1, emg sub-amp 2.
    """
    tin = next(Path(block).glob("*.tin"), None)
    if tin is None:
        return {}
    try:
        with zipfile.ZipFile(tin) as z:
            inner = next(n for n in z.namelist() if n.endswith(".synexpz"))
            with zipfile.ZipFile(io.BytesIO(z.read(inner))) as z2:
                experiment = json.loads(z2.read(next(n for n in z2.namelist() if n.endswith(".synexp"))))
    except (StopIteration, zipfile.BadZipFile, json.JSONDecodeError, KeyError):
        return {}
    gizmos = []
    def walk(o): #settings sit in nested lists/dicts, every gizmo has an ObjectState with a StateHash
        if isinstance(o, dict):
            state = o.get("ObjectState|State")
            if isinstance(state, dict):
                gizmos.append(state.get("QHash<QString,QVariant>|StateHash", {}))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(experiment)
    pz5 = next((g for g in gizmos if str(g.get("objectID", "")).startswith("PZn")), None)
    if pz5 is None:
        return {}
    out = {}
    for n in range(1, 5): #a PZ5 has up to 4 sub-amps
        send, count = pz5.get(f"sendLinkID{n}"), pz5.get(f"saNumChan{n}")
        if not pz5.get(f"saEnab{n}") or pz5.get(f"sendDisab{n}") or pz5.get(f"sendNumChan{n}") != count:
            continue
        for g in gizmos:
            if (g.get("storeID") and g.get("sourceLinkID") == f"{pz5['objectID']}.{send}" and g.get("sourceFirstChan") == 1
                    and g.get("sourceNumChan") == count and not g.get("storeDisab")):
                out[n] = (g["storeID"].replace("/", "_"), count) #same renaming tdt does ("PeA/" -> "PeA_")
    return out


def _bv_impedance(vhdr):
    #brainvision writes one "Impedance [KOhm] at hh:mm:ss (recording started at hh:mm:ss)" block in the .vhdr, then "  1 ECOGA1: 24" lines
    lines = vhdr.read_text(errors="ignore").splitlines()
    names = re.findall(r"^Ch\d+=([^,]*),", "\n".join(lines), re.M)
    out = {"channel":[], "time":[], "frequency":[], "impedance":[], "connected":[]}
    for start, line in enumerate(lines):
        clocks = re.findall(r"\d\d:\d\d:\d\d", line)
        if not line.startswith("Impedance [") or len(clocks) < 2:
            continue
        check, began = (datetime.strptime(c, "%H:%M:%S") for c in clocks[:2])
        time = (check - began).total_seconds()
        time = time + 86400 if time > 43200 else time #checks are minutes before the recording; a big positive gap means midnight passed
        for row in lines[start + 1:]:
            match = re.match(r"^\s*\d*\s*([^:]+):\s*(.*)$", row)
            if not row.strip() or match is None:
                break
            name, value = match.group(1).strip(), match.group(2).strip()
            if name in names:
                number = names.index(name) + 1
            elif name.upper().startswith("REF"): #REF_ref, REF_ECOG_A1 ... = whatever was the reference
                number = 0
            else: #GND, and a few channels that aren't in the recording
                continue
            is_number = re.fullmatch(r"-?\d+(\.\d+)?", value) is not None
            out["channel"].append(number)
            out["time"].append(time)
            out["frequency"].append(np.nan)
            out["impedance"].append(float(value) if is_number else np.nan)
            out["connected"].append(int(is_number)) #"Out of Range!" and "Disconnected!" are the only words seen
    return {"dbs":_as_arrays(out)} if out["channel"] else {}


def _as_arrays(rows):
    #plain numpy arrays so it saves to .h5 and opens in MATLAB
    return {"channel":np.asarray(rows["channel"], dtype=np.int64), "time":np.asarray(rows["time"], dtype=np.float64),
            "frequency":np.asarray(rows["frequency"], dtype=np.float64), "impedance":np.asarray(rows["impedance"], dtype=np.float64),
            "connected":np.asarray(rows["connected"], dtype=np.int8)}

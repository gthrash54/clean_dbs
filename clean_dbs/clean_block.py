import warnings
from pathlib import Path
from clean_dbs.read_tdt import tdt_windows, WINDOW_MB
from clean_dbs.read_bv import bv_windows
from clean_dbs.save_h5 import open_h5, StreamWriter, write_rest, stale_reason
from clean_dbs.block_info import output_path, describe_block, recording_condition, source_path
from clean_dbs.flags import make_flags, Saturation
from clean_dbs.datatypes import stream_type, event_type
from clean_dbs.channel_map import channel_table, connected_contacts, leads_implanted, leads_connected
import pandas as pd
from clean_dbs.stim_events import stim_events, bv_clock_time
from clean_dbs.markers import label_bv
from clean_dbs.coords import contact_coords, coord_status, RECON_COHORTS
from clean_dbs.dnu import dnu_channels
from clean_dbs.cohorts import cohort_of
from clean_dbs import privacy

def clean_block(block_dir, out_root, window_mb=WINDOW_MB, force=False):
    """
    Inputs:
    block_dir is the block to process: a tdt block folder (from find_tdt_blocks) or a .vhdr file (from find_bv_recordings)
    out_root is the top of the cleaned set where outputs go
    window_mb is about how many MB of samples to hold at once; the recording is read, checked and written one window at a time,
      so memory stays around 0.5 GB whatever the file size (a 13 GB brainvision file included), and the .h5 is the same for any window_mb
    force: True cleans it again even if a current .h5 is already there
    Returns: path of the saved file, out_root/<patient>/<stage>/<block>.h5 (uh/UG) or out_root/<patient>/<session>/<block>.h5 (other cohorts)
    """
    block_dir = Path(block_dir)
    out_root = Path(out_root)
    privacy.use_out_folder(out_root) #names in <out_root>/_cases/renames.csv are written as their new names everywhere below
    out_path = output_path(block_dir, out_root) #works out patient/stage/name, see block_info.py
    if out_path.exists() and not force and not stale_reason(out_path): #done by this CLEAN_VERSION; broken or older files are cleaned again
        return out_path
    if block_dir.suffix.lower() == ".vhdr": #brainvision is a .vhdr file, tdt is a block folder
        r, windows = bv_windows(block_dir, window_mb)
    else:
        r, windows = tdt_windows(block_dir, window_mb)
    case_targets = _case_targets(out_root, block_dir) #lead 1 / lead 2 targets for this patient from the OR sheets, and that row's status
    case_row = _case_row(out_root, block_dir)
    implanted = leads_implanted(case_row) #a lead that never went in is labeled lead<N>_not_implanted
    connected = leads_connected(case_row) #one that went in but was never plugged in is labeled lead<N>_not_connected
    about = describe_block(block_dir)
    #recon, no_recon or no_contacts_csv, kept in the file so a missing contacts.csv can't pass for a patient with no recon
    coords_found = coord_status(out_root, about["patient"])
    if coords_found == "no_contacts_csv" and cohort_of(about["patient"]) in RECON_COHORTS:
        _warn_no_contacts(out_root)
    is_bv = block_dir.suffix.lower() == ".vhdr"
    for name, s in r["streams"].items(): #tag each stream lfp / mer / emg / audio / stim_monitor, see datatypes.py (read_bv already typed its streams)
        if "type" not in s:
            s["type"], s["subtype"] = stream_type(s["device"])
        table = channel_table(block_dir, name, s["device"], len(s["channels"]), *case_targets[:3], implanted=implanted, connected=connected) #what each channel is (see channel_map.py)
        if table is not None:
            s["channel_names"] = table.pop("name")
            for key, values in table.items():
                s[f"channel_{key}"] = values
            #where each lead contact was (Brainlab ACPC mm) during this recording, from _cases/contacts.csv; NaN = not known or not in yet
            about_block = describe_block(block_dir)
            xs, ys, zs, positions, space = contact_coords(out_root, about_block["patient"], about_block["block"],
                                                          s["channel_group"], s["channel_contact"], about_block.get("stage", ""))
            if space:
                s["channel_x"], s["channel_y"], s["channel_z"], s["channel_position"] = xs, ys, zs, positions
                s["coord_space"] = space
        #channels the lab marked do not use in this recording (bad contacts, contacts still in the cannula), from _cases/dnu_channels.csv;
        #dnu_list says whether the recording was on the list at all ("none" = not listed, so nothing is known), see dnu.py
        listed_as = describe_block(block_dir)
        s["channel_dnu"], listed = dnu_channels(out_root, listed_as["patient"], listed_as.get("stage", ""), listed_as["block"],
                                                s.get("channel_names"), s["channels"], s["type"])
        s["dnu_list"] = listed or "none"
        if s["type"] == "mer" and not is_bv: #AlphaOmega MER comes into the TDT ADC after the Neuro Omega analog out (total gain 4000, manual
            s["analog_gain"] = 4000.0        #v1.4.6 p183/202): electrode volts = data / analog_gain. Data left as recorded; routing unconfirmed
    #the video stays outside the .h5 (too big, and video tools cant read it from inside); record where it is, Cam1 lines frames up with the data
    videos = sorted(block_dir.glob("*.avi")) if block_dir.is_dir() else []
    r["info"]["video"] = ";".join(source_path(block_dir) + "/" + privacy.rename_file(v.name) for v in videos) #from the study ID down, like source
    #what the recording is, from the folder and file names (see cohorts.py), so cohorts can be pooled by condition later; blank = names don't say
    r["info"]["condition"], r["info"]["state"] = recording_condition(block_dir)
    r["info"]["case_status"] = case_targets[3] #how far to trust the targets above: match, conflict, sheet_only, unchecked, missing (no row)
    r["info"]["coord_status"] = coords_found
    if about["block"] != (block_dir.stem if is_bv else block_dir.name): #renamed in renames.csv; the reader took it from the file name
        r["info"]["block"] = about["block"]
    for key in ("cohort", "session", "parent_id"): #only other cohorts have these
        if key in about:
            r["info"][key] = about[key]
    r["info"].setdefault("cohort", cohort_of(about["patient"]) or "") #uh/UG too, lowercase like cohorts.py (tdt doesn't say)
    saturation = {name:Saturation(s) for name, s in r["streams"].items()}
    with open_h5(out_path) as f: #written as <block>.h5.part, renamed to .h5 only once everything below is in
        writer = StreamWriter(f, r["streams"])
        for window in windows: #one stretch of every stream at a time, never the whole recording
            writer.write(window)
            for name, (_, x) in window.items():
                saturation[name].add(x)
            window = x = None #so this window is freed before the next one is read
        writer.close()
        #saturation, new segments, dnu, stored as notes, the data isn't changed
        r["flags"] = make_flags(r, block_dir, {name:s.result() for name, s in saturation.items()})
        if is_bv: #marker codes -> stim / task / marker with their meaning (markers.py), stim table, and the stimulator files next to the recording
            labeled = label_bv(r, block_dir)
            r["stim"], r["sidecars"] = labeled["stim"], labeled["sidecars"]
            privacy.scrub_result(r) #sidecar file names (in /sidecars, the stim table and marker evidence) as renames.csv says
        else:
            for name, e in r["epocs"].items():
                e["type"] = event_type(name)
            for name, sc in r["scalars"].items():
                sc["type"] = event_type(name)
            r["stim"] = stim_events(r, block_dir) #one row per stim pulse: time, voice, channel, amplitude, period, count
        if "new_segments" in r["flags"]: #brainvision: the same pulse times on the real clock, with the missing time at real gaps added back
            r["stim"]["time_clock"] = bv_clock_time(r["stim"]["time"], r["flags"]["new_segments"])
        for name, check in connected_contacts(block_dir).items(): #impedance checks from the rig, connected=0 means open/disconnected
            r["flags"][f"impedance_{name}"] = check
        write_rest(f, r, block_dir)
    return out_path


#out folders already warned about, so a run over hundreds of uh blocks says it once, not once per block
_warned_no_contacts = set()


def _warn_no_contacts(out_root):
    #uh patients have recons, so cleaning them without _cases/contacts.csv quietly leaves every file with no coordinates
    if str(out_root) in _warned_no_contacts:
        return
    _warned_no_contacts.add(str(out_root))
    warnings.warn("no _cases/contacts.csv in the out folder, so the uh files from this run get no contact coordinates "
                  "(info coord_status = no_contacts_csv). Copy contacts.csv into _cases and clean them again with --force", stacklevel=3)


def _case_targets(out_root, block_dir):
    #lead 1 / lead 2 targets (e.g. STN, GPi), the ecog entry and the cross-check status from <cleaned root>/_cases/cases.csv
    #blank targets and "missing" if the case isn't in the table (clean.py stops before this when cases.csv isn't there, unless --no-cases)
    row = _case_row(out_root, block_dir)
    if row is None:
        return ("", "", "", "missing")
    status = row.get("status", "").strip() #status from scripts/crosscheck_cases.py, blank = never checked
    return (row["target1"], row["target2"], row["ecog"], status or "unchecked")


def _case_row(out_root, block_dir):
    #this patient's row of <cleaned root>/_cases/cases.csv as a dictionary, None if there's no table or no row
    cases = Path(out_root) / "_cases" / "cases.csv"
    if not cases.is_file():
        return None
    table = pd.read_csv(cases, dtype=str, keep_default_na=False)
    #lowercase both sides, Box has UG3003 where the tables have ug3003; patient, not parent.parent, since other cohorts sit deeper
    row = table[table["study_id"].str.strip().str.lower() == describe_block(block_dir)["patient"].lower()]
    return row.iloc[0].to_dict() if len(row) else None

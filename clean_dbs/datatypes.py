from pathlib import Path
import h5py

#device (Synapse gizmo name, or "brainvision") -> (type, subtype)
#built from every device in the inventory on 2026-10-11; anything not listed is tagged "unknown"
DEVICE_TYPES = {
    "ecogdbs":("lfp", "ecog+lead"), #strip and lead contacts in one store, splitting them needs the OR data sheets
    "ecog":("lfp", "ecog"),
    "dbs1":("lfp", "lead"),
    "brainvision":("lfp", ""), #stage 2 is the lead; some ug stage 1 brainvision files are unclear, so no subtype yet
    "raw":("mer", "raw"),
    "spk":("mer", "spike"),
    "spike":("mer", "spike"),
    "macraw":("mer", "macro_raw"),
    "microlfp":("mer", "micro_lfp"),
    "macrolfp":("mer", "macro_lfp"),
    "maclfp":("mer", "macro_lfp"),
    "emg":("emg", ""),
    "microphone":("audio", ""),
    "MonA":("stim_monitor", ""),
    "RSn1":("unknown", ""), #RS4 data streamer, 16 ch in uV, only on uh3003, nobody knows what it is yet
}


def stream_type(device):
    """
    Input: device name (the "device" a reader puts on each stream)
    Returns: (type, subtype), e.g. ("mer", "spike"); ("unknown", "") if the device isn't in DEVICE_TYPES
    """
    return DEVICE_TYPES.get(device, ("unknown", ""))


def load_type(h5_path, data_type, subtype=None):
    """
    Loads just the streams of one type from a cleaned .h5 file, e.g. all the mer streams.
    Input:
    h5_path: a cleaned file from clean.py
    data_type: "lfp", "mer", "emg", "audio", "stim_monitor" or "unknown"
    subtype: optional, e.g. "spike" to get only the spike-band mer streams
    Returns: {stream name: {"data", "fs", "start_time", "channels", "device", "type", "subtype", "units"}}
    """
    out = {}
    with h5py.File(Path(h5_path), "r") as f:
        for name, g in f["streams"].items():
            if g.attrs.get("type") == data_type and (subtype is None or g.attrs.get("subtype") == subtype):
                out[name] = {"data":g["data"][:], "channels":g["channels"][:], **{key:g.attrs[key] for key in g.attrs}}
    return out

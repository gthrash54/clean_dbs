import re
from pathlib import Path
import h5py

#device (Synapse gizmo name, or "brainvision") -> (type, subtype)
#built from every device in the inventory on 2026-10-11; anything not listed is tagged "unknown"
DEVICE_TYPES = {
    "ecogdbs":("lfp", "ecog+lead"), #strip and lead contacts in one store, splitting them needs the OR data sheets
    "ecog":("lfp", "ecog"),
    "dbs1":("lfp", "lead"),
    "brainvision":("lfp", ""), #only for read_bv(split=False); the normal read_bv types each stream from its channel names (bv_channels.py)
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
    "RSn1":("unknown", ""), #RS4 data streamer, 16 ch in uV, only on uh3002 (all 14 local blocks), nobody knows what it is yet
}


#event stores (tdt epocs/scalars, brainvision markers) -> type
#stim = one event per pulse; stim_params = a stimulator setting logged next to the pulses (count, amplitude, delay, channel, ...), not a pulse,
#so index.n_stim_events and the /stim table only count the pulse stores. Store names seen on the 408 local blocks (2026-10-12):
#PeA_/PeB_/PeC_ (+ PerB on uh3001) and the settings AmpA-C, CntA-C, CnA_, CnB_, DelA-C, ChnA-B
EVENT_TYPES = {
    "Cam1":"video", #one event per camera frame, data = frame number, lines frames up with the neural clock
    "PeA_":"stim", "PeB_":"stim", "PeC_":"stim", "PeD_":"stim", #PeA_ = one event per pulse on voice A, its data is the period in ms
    "PerB":"stim", #uh3001 names voice B's pulse store PerB (crosstalk1/2): same onsets and period values as PeB_ elsewhere
    "eS1p":"stim_params", #24 rows = period, count, amp (uA), duration (ms), delay (ms), channel for voices A-D, from the lab's Synapse experiment; period row lags one train, use PeA_ data instead
    "Note":"other", "Tick":"other",
}
#any other eStim setting saved as its own store: Amp/Am, Cnt/Cn, Del/De, Dur/Du, Chn/Ch, Per + voice letter, with or without the "_" tdt adds for "/"
SETTING_RE = re.compile(r"(Amp|Am|Cnt|Cn|Del|De|Dur|Du|Chn|Ch|Per)[A-D]_?")


def event_type(name, is_bv=False):
    """
    Input: event store name (e.g. Cam1, PeA_, R15), is_bv for brainvision markers
    Returns: video, stim, stim_params, marker or other
    brainvision markers (R15, S13, ...) mean different things per cohort: clean_block types them with markers.label_bv (stim / task / marker + meaning);
    this fallback calls them all "marker"
    """
    if is_bv:
        return "other" if name == "NewSegment" else "marker"
    if name in EVENT_TYPES:
        return EVENT_TYPES[name]
    return "stim_params" if SETTING_RE.fullmatch(name) else "other"


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
    data_type: "lfp", "mer", "emg", "audio", "stim_monitor", "eeg", "accel", "joystick", "force", "aux", "unused" or "unknown"
    subtype: optional, e.g. "spike" to get only the spike-band mer streams
    Returns: {stream name: {"data", "fs", "start_time", "channels", "device", "type", "subtype", "units", "rail", ...}} plus, when the
             stream has them, what each channel is, like clean_dbs_load.m gives: channel_names, channel_group, channel_target,
             channel_contact, channel_source, ... as lists of text, and the contact coordinates channel_x, channel_y, channel_z as
             float arrays in mm (NaN = unknown) with channel_position and the stream's coord_space
    """
    out = {}
    with h5py.File(Path(h5_path), "r") as f:
        for name, g in f["streams"].items():
            if g.attrs.get("type") == data_type and (subtype is None or g.attrs.get("subtype") == subtype):
                out[name] = {"data":g["data"][:], "channels":g["channels"][:], **{key:g.attrs[key] for key in g.attrs}}
                for key in g: #per-channel labels and coordinates, saved as one dataset each next to data (see save_h5)
                    if key.startswith("channel_"):
                        is_text = h5py.check_string_dtype(g[key].dtype) is not None
                        out[name][key] = g[key].asstr()[:].tolist() if is_text else g[key][:]
    return out

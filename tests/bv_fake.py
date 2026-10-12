from pathlib import Path
import numpy as np


def write_bv(folder, name, channels, data, fs, markers, dated=True):
    """
    Writes a tiny made-up brainvision recording (.vhdr, .eeg float32 multiplexed, .vmrk) that MNE can read.
    Input: folder, name (file stem), channels [(name, unit)], data (channels x samples, in the units written), fs (Hz),
           markers [(type, description, sample 1-based)], dated (False leaves the date off the New Segment, like a recording without one)
    Returns: path of the .vhdr
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    lines = ["Brain Vision Data Exchange Header File Version 1.0", "[Common Infos]", "Codepage=UTF-8", f"DataFile={name}.eeg",
             f"MarkerFile={name}.vmrk", "DataFormat=BINARY", "DataOrientation=MULTIPLEXED", f"NumberOfChannels={len(channels)}",
             f"SamplingInterval={1e6 / fs:g}", "[Binary Infos]", "BinaryFormat=IEEE_FLOAT_32", "[Channel Infos]"]
    lines += [f"Ch{i + 1}={n},,1,{u}" for i, (n, u) in enumerate(channels)]
    (folder / f"{name}.vhdr").write_text("\n".join(lines) + "\n", encoding="utf-8")
    np.asarray(data, dtype="<f4").T.tofile(folder / f"{name}.eeg")
    stamp = ",20240101120000000000" if dated else ""
    mk = [f"Mk1=New Segment,,1,1,0{stamp}"] + [f"Mk{i + 2}={t},{d},{s},1,0" for i, (t, d, s) in enumerate(markers)]
    (folder / f"{name}.vmrk").write_text("Brain Vision Data Exchange Marker File, Version 1.0\n[Common Infos]\nCodepage=UTF-8\n"
                                         f"DataFile={name}.eeg\n[Marker Infos]\n" + "\n".join(mk) + "\n", encoding="utf-8")
    return folder / f"{name}.vhdr"


def pulse_markers(starts, width_samples, start="S 15", end="S 13"):
    #one TTL per start sample: the start code, then the end code width_samples later
    out = []
    for s in starts:
        out += [("Stimulus", start, int(s)), ("Stimulus", end, int(s + width_samples))]
    return out


def write_mcs(path, rows_ch1, rows_sync, mode="current", sync_channel=9):
    """
    Writes a made-up MC_Stimulus II ASCII program (format 3): channel 1 = stim rows, sync_channel = TTL rows.
    Input: path, rows_ch1 / rows_sync = lists of (value, time_us, value, time_us, value, time_us, repeat)
    Returns: path
    """
    out = ["Multi Channel Systems MC_Stimulus II", "ASCII import Version 1.10", "", "channels:\t8", f"output mode:\t{mode}", "format:\t3", ""]
    for ch in range(1, sync_channel + 1):
        rows = rows_ch1 if ch == 1 else rows_sync if ch == sync_channel else []
        out += [f"channel:\t{ch}", "value\ttime\tvalue\ttime\tvalue\ttime\trepeat"] + ["\t".join(str(int(v)) for v in r) for r in rows] + [""]
    Path(path).write_text("\r\n".join(out))
    return Path(path)

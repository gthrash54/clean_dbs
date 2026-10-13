import re
from pathlib import Path
import numpy as np
import pandas as pd
from clean_dbs.cohorts import ID_IN_NAME_RE

#Some folder and file names on Box carry a person's first name (a recording saved as the patient's name, a subfolder named after
#whoever ran it). Those names would end up in the cleaned folders, file names, the "source" and "block" info in every file, and the
#reports. So every name that goes into an output has to be made only of the words below (what was recorded, how, when), or of words a
#person added to <out>/_cases/allowed_words.txt, or be renamed in <out>/_cases/renames.csv. make_manifest.py refuses to go on otherwise
#and lists what it couldn't clear in <out>/_cases/names_to_check.csv. When in doubt it refuses: a word it doesn't know is never let through.
#study IDs (DLPFC045, uh3021) are taken out before the check, numbers and single letters are fine. Add words here only if they describe
#a recording; never a name, and never a word that is mostly used as a first name.
WORDS = frozenset("""
uh ug dlpfc dlpf bi pp mc bv orbv bpp bb chr pt sep pns tm tmscs lfp ecog burst dys sensory pr motif
sdc or stage surgery surgical sedated not ga awake asleep anesth anesthesia anesthetized precedex propofol intraop preop postop op pre post dbs stn gpi
vim zi snr thal thalamus cortex cortical ctx motor premotor frontal parietal lfp lfps ecog eeg ieeg hdeeg emg ecg ekg eog mer mers
hd scalp strip lead leads contact contacts comtact electrode electrodes macro micro macrostim microstim stim stims stimulation stimulus stimuli stimulator nostim teststim corticalstim prestim
poststim stimfile dbsstim merstim reversalstim erpstim erpassaystim wbdbs wdbs pnsdbs rest resting baseline baselinerest baselinetask baselinelfp eyes eyesopen eyesclosed resteyesopen resteyesclosed open closed cross
wcross fixation task tasks tssk practice prac axcpt simon simonv simontask motif joystick reach choice go nogo gonogo flanker stroop nback oddball itbs ctbs
tbs theta burst bursts burstpp burstspp burstsandpp and pairedpulse pairedpulses paired pulse pulses pair pairs single singlepulse double train trains ramp freqramp sdcurve sd
curve erp erps erpassay assay mep threshold thresh threhold mepthreshold probing polarity reversal rev arie ariebest arieworst best worst home settings homesettings setting pseudo
psuedo pseudomonopolar homepseudomonopolar monopolar bipolar multipolar mono therapeutic ther window capsular capsularmonopolar capsularpseudomonopolar bestmonopolar worstmonopolar bestpseudomonopolar worstpseudomonopolar directional directionality dir initial peripheral passive active
movement movements move voluntary lesion lesionlfp buzz buzzy bee hand oc off on meds med offmeds onmeds predbs postdbs postmer hz khz ma ua
mv uv us ms sec min hr day days month months mo hfs lf high low freq frequency highfreq lowfreq amp amplitude width pw
isi strength duration range vary volts fixed static dynamic left right lt rt bilateral side hand arm forearm elbow leg mouth face foot wrist
finger bicep tremor data raw folder contents of copy new old run trial trials block blck session part test testing clinicaltesting misc results ttl
ppttl protocol stimuli ordered randomized voltages impedance impedence impedencetest impedencecheck impedencetab ecogimpedance check clinical increment cam microphone crosstalk code date time table survey metronome
rainbow passage cognitive cog implicit collisions induction during after before from by only no first second repeat rep redo real actual modified processed annotated
reposition shift search found find everything force cart bottom stem target loc place me co cat er ps pn db pi sp stg inst
ecogonly microrecord cont exp laptop nochronos vhdr max multi state
brain aim traj trajectory depth events channels struct ring record files final neuro omega alpha ao medtronic incrementdbs
good bad instdelreach neuromega multifreq multipw afterdbs notask nofixation restingstateprestimulation restingstatestimulation restwcross pairpulse trolling
""".split())

CASES = "_cases"
RENAMES = "renames.csv"              #study_id, original, renamed: what to write instead of a folder or file name (renamed = original keeps it)
ALLOWED = "allowed_words.txt"        #one word per line, words that describe recordings but aren't in WORDS
TO_CHECK = "names_to_check.csv"      #written by make_manifest.py: everything it couldn't clear, for a person to fill in
WITHHELD = "<withheld>"              #what a name that isn't cleared turns into when printed or reported

_active = {"key":None, "renames":{}, "allowed":frozenset()}


def use_out_folder(out_root):
    """
    Loads the rename map and the allowed words from <out_root>/_cases; describe_block, source_path and output_path use them from then on.
    Input: out_root, the cleaned set's top folder (None or a folder without the files = no renames, built-in words only)
    Returns: (number of renames, number of allowed words)
    Raises ValueError when renames.csv has a renamed value that doesn't pass the check itself (a typo could put a name back in)
    """
    folder = Path(out_root) / CASES if out_root is not None else None
    files = [folder / RENAMES, folder / ALLOWED] if folder else []
    key = tuple((str(f), f.stat().st_mtime_ns) for f in files if f.is_file())
    if key == _active["key"]:
        return len(_active["renames"]), len(_active["allowed"])
    allowed = set()
    if folder and (folder / ALLOWED).is_file():
        for line in (folder / ALLOWED).read_text(errors="ignore").splitlines():
            #each run of letters is one allowed word as written, never split ("InstDelReach" allows instdelreach, not "del" on its own)
            allowed.update(run.lower() for run in _runs(line.split("#", 1)[0]))
    renames = {}
    if folder and (folder / RENAMES).is_file():
        table = pd.read_csv(folder / RENAMES, dtype=str, keep_default_na=False)
        missing = {"original", "renamed"} - set(table.columns)
        if missing:
            raise ValueError(f"{folder / RENAMES} needs the columns study_id, original, renamed (missing {sorted(missing)})")
        for row in table.to_dict("records"):
            original, renamed = row["original"].strip(), row["renamed"].strip()
            if not original or not renamed: #a row nobody filled in yet doesn't count
                continue
            if re.search(r"[/\\]", renamed) or renamed.strip(". ") == "":
                raise ValueError(f"{RENAMES}: a renamed value has a / or \\ in it or is only dots (row for study {row.get('study_id', '')!r})")
            if unknown_words(renamed, allowed):
                raise ValueError(f"{RENAMES}: a renamed value doesn't pass the name check itself (row for study {row.get('study_id', '')!r}); "
                                 "use only words that describe the recording")
            renames[(row.get("study_id", "").strip().upper(), original.lower())] = renamed
    _active.update(key=key, renames=renames, allowed=frozenset(allowed))
    return len(renames), len(allowed)


def reset():
    #back to no renames and built-in words only (tests, or before switching to another out folder)
    _active.update(key=None, renames={}, allowed=frozenset())


def rename(part, study_id=""):
    """
    Input: part, one folder name or file stem as on Box; study_id, whose folder it's in ("" = any)
    Returns: what renames.csv says to write instead, or part unchanged; a row for this study ID wins over a row with no study ID
    """
    renames = _active["renames"]
    if not renames:
        return part
    key = str(part).strip().lower()
    return renames.get((str(study_id).upper(), key), renames.get(("", key), part))


def rename_file(name, study_id=""):
    #a file name: the stem goes through rename, the extension stays
    stem, dot, ext = str(name).rpartition(".")
    if not dot or not stem:
        return rename(name, study_id)
    return rename(stem, study_id) + "." + ext


def words_of(text):
    """
    Input: a folder name or file stem
    Returns: its words in lowercase, without study IDs, numbers and camera numbers ("HandOCAndBuzz2" -> hand, oc, and, buzz;
        "DLPFC045_rest" -> rest; "lfp1_Cam1" -> lfp)
    """
    return [w for run in _runs(text) for w in _splits(run)[1]]


def _runs(text):
    #runs of letters, after taking out study IDs and tdt's Cam1/Cam2 video names
    text = re.sub(r"(?i)(?<![a-z])cam(?=\d)", " ", ID_IN_NAME_RE.sub(" ", str(text)))
    return re.findall(r"[^\W\d_]+", text)


def _splits(run):
    #ways to read one run of letters as words: as it is, split the camel case way (HandOCAndBuzz -> Hand OC And Buzz), split at every
    #change of case (BRAINaim -> BRAIN aim). Letters outside a-z stay one word, so an accented name is never split into known bits
    if not run.isascii():
        return [[run.lower()]] * 2
    camel = re.findall(r"[A-Z]+(?![a-z])|[A-Z]?[a-z]+", run)
    case = re.findall(r"[A-Z]+|[a-z]+", run)
    return [[w.lower() for w in way] for way in ([run], camel, case)]


def unknown_words(text, allowed=None):
    """
    Input: text, a folder name or file stem; allowed, extra words (defaults to the ones loaded by use_out_folder)
    Returns: the words in it that aren't known (empty list = nothing in it can be a name). A run of letters passes when one way of
        reading it (see _splits) is all known words; single letters always pass
    """
    allowed = _active["allowed"] if allowed is None else allowed
    out = []
    for run in _runs(text):
        ways = [[w for w in way if len(w) > 1 and w not in WORDS and w not in allowed] for way in _splits(run)]
        if all(ways):
            out += ways[1]
    return out


def is_clear(text):
    return not unknown_words(text)


def mask(path_text):
    """
    Input: a path from the study ID down ("DLPFC045/Misc/x.vhdr"), for printing or a report
    Returns: the same path with every folder or file name that doesn't pass the check replaced by <withheld>
    """
    parts = str(path_text).replace("\\", "/").split("/")
    return "/".join(p if is_clear(Path(p).stem if i == len(parts) - 1 else p) else WITHHELD for i, p in enumerate(parts))


def scrub(text):
    """
    Input: free text that can hold raw file or folder names (an error message)
    Returns: the text with every name in renames.csv swapped for its new name
    """
    #longest first, so a file named after a folder ("<folder>_results") is swapped whole, not just its folder part
    for (_, original), renamed in sorted(_active["renames"].items(), key=lambda item: -len(item[0][1])):
        text = re.sub(rf"(?<![A-Za-z]){re.escape(original)}(?![A-Za-z])", lambda m: renamed, str(text), flags=re.IGNORECASE)
    return text


def scrub_result(result):
    """
    Input: the dictionary clean_block writes (reader output plus labels), changed in place
    Swaps names from renames.csv in the text that came from file names: sidecar names and keys, the stim table's sidecar column, and the
    evidence notes on the markers (the file contents themselves aren't touched)
    """
    if not _active["renames"]:
        return result
    for e in result.get("epocs", {}).values():
        for key, value in e.items():
            if isinstance(value, str):
                e[key] = scrub(value)
    if "sidecar" in result.get("stim", {}):
        result["stim"]["sidecar"] = np.asarray([scrub(v) for v in result["stim"]["sidecar"]], dtype=object)
    sidecars = {}
    for key, entry in result.get("sidecars", {}).items():
        renamed = scrub(entry["file"])
        if renamed != entry["file"]: #the group name is made from the file name the way sidecars._key does it ("x_1.1.txt" -> x_1_1)
            entry["file"] = renamed
            key = re.sub(r"[^A-Za-z0-9_]+", "_", Path(renamed).stem).strip("_") or "sidecar"
        while key in sidecars:
            key += "_"
        sidecars[key] = entry
    if "sidecars" in result:
        result["sidecars"] = sidecars
    return result

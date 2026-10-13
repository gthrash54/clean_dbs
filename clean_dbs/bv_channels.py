import re
from pathlib import Path
from clean_dbs.cohorts import cohort_of as study_cohort

#what each brainvision channel is, from the name the recorder wrote in the .vhdr. Built from every .vhdr on Box (1,940 headers,
#24 cohorts, survey 2026-10-12). Names are matched after lowercasing and turning spaces into "_" ("EMG FDI_Ipsi" -> "emg_fdi_ipsi").
#Each pattern -> (type, subtype). Order matters, the first match wins. Anything not matched is "unknown": we never guess.
#types: lfp (subtype lead / ecog), eeg (scalp), emg, audio, mer, accel, joystick, force, aux, unused, unknown
MUSCLES = {"fdi":"first dorsal interosseous", "fcr":"flexor carpi radialis", "fcu":"flexor carpi ulnaris", "ecr":"extensor carpi radialis",
           "ecu":"extensor carpi ulnaris", "oo":"orbicularis oris", "orbicularis":"orbicularis", "gastroc":"gastrocnemius", "leggastroc":"gastrocnemius",
           "ta":"tibialis anterior", "tib_ant":"tibialis anterior", "legquad":"quadriceps", "bicep":"biceps", "armbicep":"biceps", "tricep":"triceps",
           "thenar":"thenar", "armext":"arm extensors", "armflex":"arm flexors", "ext":"extensors", "flex":"flexors", "emg":""}
SIDES = {"contra":"contra", "conta":"contra", "ipsi":"ipsi", "left":"left", "right":"right"} #"conta" is a typo in some headers
MUSCLE = "(?P<muscle>" + "|".join(sorted(MUSCLES, key=len, reverse=True)) + ")"
SIDE = "(?P<side>contra|conta|ipsi|left|right)"
PATTERNS = [
    (r"", ("unused", "empty_slot")), #blank name = nothing plugged in (one BPP header)
    (r"no_(micro|macro|spike|raw|lfp)_?(raw|spk|spike)?", ("unused", "empty_mer_slot")), #NO_MICRO_RAW, NO_LFP ...: MER inputs labeled not connected
    (r"(dbs|lfp|db_s)_?(?P<contact>\d{1,2})(_(?P<hemi>[lr]))?|lfp", ("lfp", "lead")), #DBS1-16, LFP1-8, DBS1_L ("DB S8" is a typo for DBS8)
    (r"ecog_?(?P<a>[ap])?_?(?P<contact>\d{1,2})_?(?P<b>[ap])?", ("lfp", "ecog")), #ECOG1, ECOG_A1, ecog_a1, ecog_1a, ECOGp1
    (r"(fp|af|f|fc|ft|c|cp|t|tp|p|po|o|i)(z|\d{1,2})", ("eeg", "scalp")), #10-10 electrode names (Cz, FCz, AF7, Iz, "Y1 AF7" ...)
    (r"scalpeeg|(left|right)_scalp", ("eeg", "scalp_unlabeled")), #scalp electrode, position not written
    (r"mic|microphone", ("audio", "mic")),
    (r"(ao_)?micro_?(raw|spike|spk)?_?\d*|ao_raw\d*_micro|microraw", ("mer", "micro")), #AO = Alpha Omega, the MER system
    (r"(ao_)?macro_?(raw|spike|spk|lfp)?_?\d*|macroraw|macrolfp", ("mer", "macro")),
    (r"raw_?\d*", ("mer", "raw")), #RAW_01 next to Macro_LFP01 and SPK01: the three Alpha Omega outputs
    (r"spk_?\d*|spikes?", ("mer", "spike")),
    (r"mer|ao_\d+", ("mer", "")), #MER, ao_1: a MER system output, band not named
    (r"acc(elerometer)?_?(?P<axis>[xyz])", ("accel", "")),
    (r"joy_?(?P<axis>[xy])", ("joystick", "")),
    (r"dynamometer|force|grip", ("force", "")),
    (r"photo_?sensor", ("aux", "photo_sensor")), #light sensor on the task screen
    (r"aux\d+", ("aux", "")), #actiCHamp AUX inputs, nothing in the header says what was plugged in
    (rf"(emg_)?{MUSCLE}\d?(_?(?P<part>ext|flex|lower))?_?{SIDE}?\d*", ("emg", "")),
    (rf"{SIDE}_{MUSCLE}", ("emg", "")), #Contra_FCR
]
#names seen on Box that stay unknown on purpose, so nobody "fixes" them by guessing:
#  Ch1..Ch64 and bare numbers 1..32 = the amplifier input number, renamed by nobody (BI, PT, TM, PNS, DYS, DLPFC, MC, CHR, LFP, ECOG, PP, TMSCS)
#  G1..G32 = actiCAP holder numbers (DLPFC061 preop); Arm, Leg, Mouth, Foot, knee_right = a body site, could be EMG, a sensor or a trigger;
#  SCS1-8 (TMSCS); EP1-8 (PP); FRONTAL1/2, Ch33Frontal (frontal what?); FCI (not a muscle abbreviation we know)
UNKNOWN_NOTES = {r"(ch_?)?\d+":"unlabeled amplifier input", r"g\d+":"actiCAP holder number, position not in header",
                 r"(arm|leg|mouth|foot|knee)(_\w+)?":"body site only", r"scs\d+":"SCS contact (not documented)", r"ep\d+":"not documented"}

#per-cohort mappings for channels the header leaves generic, only where something besides the name says what they are.
#cohort -> list of (name pattern, (type, subtype), evidence). Add rows here when the lab confirms a montage; leave the rest unknown.
COHORT_OVERRIDES = {
    "BV":[(r"ch\d+", ("eeg", "scalp_unlabeled"), "BV cohort = 64-channel HD-EEG on actiCHamp (Ch1-64 + Aux); cap position of each ChN isn't in the header, needs the lab's cap layout")],
    "ORBV":[(r"ch\d+", ("eeg", "scalp_unlabeled"), "ORBV cohort = 64-channel HD-EEG on actiCHamp (Ch1-64 + Aux); cap position of each ChN isn't in the header, needs the lab's cap layout")],
}

#stream name for each type, so a recording becomes streams like lead, ecog, eeg, emg, audio, mer, ...
STREAM_NAMES = {("lfp", "lead"):"lead", ("lfp", "ecog"):"ecog", "eeg":"eeg", "emg":"emg", "audio":"audio", "mer":"mer", "accel":"accel",
                "joystick":"joystick", "force":"force", "aux":"aux", "unused":"unused", "unknown":"unknown"}
VOLT_UNITS = {"V":1.0, "mV":1e-3, "µV":1e-6, "uV":1e-6, "μV":1e-6, "nV":1e-9} #units MNE converts to volts; anything else stays as written


def read_vhdr_channels(vhdr):
    """
    Input: a brainvision .vhdr
    Returns: list of (name, units) per channel, in file order (channel number = position + 1); units "µV" when the header leaves it blank
    """
    text = Path(vhdr).read_text(errors="ignore")
    out = []
    for m in re.finditer(r"^Ch(\d+)=([^\r\n]*)$", text, re.M):
        fields = m.group(2).split(",")
        name = fields[0].replace("\\1", ",").strip()
        units = fields[3].strip() if len(fields) > 3 and fields[3].strip() else "µV" #brainvision's default unit
        out.append((name, units))
    return out


def cohort_of(path):
    """
    Input: any path inside a study folder (e.g. .../DLPFC064/baseline.vhdr, .../uh3021/stage 2 surgery/STN_dbs1.vhdr)
    Returns: the cohort in capitals from the nearest folder above it that is a study ID (cohorts.cohort_of, so both agree on what an ID is), "" if none
    The path version of cohorts.cohort_of, which takes one folder name and gives lowercase; markers.RULES and COHORT_OVERRIDES use capitals.
    """
    for part in reversed(Path(path).parent.parts):
        cohort = study_cohort(part) #a study ID has digits; a plain "LFP" or "Stim" subfolder, DIR002 or Run-001 isn't one
        if cohort:
            return cohort.upper()
        if part.strip().upper() in ("MOTIF", "PR", "BURSTPAIRPULSES"): #the few study folders with no number
            return part.strip().upper()
    return ""


def normalize(name):
    """
    Input: a channel name as written in the .vhdr
    Returns: lowercase, spaces -> "_", with label decorations removed: a trailing " Ch7" input number ("DBS1 Ch7"), a leading actiCAP
             holder ("Y1 AF7" -> "af7"), a trailing cohort tag ("FCR_contra_BB")
    """
    n = re.sub(r"\s+", "_", name.strip().lower())
    n = re.sub(r"^(.+?)_ch_?\d+$", r"\1", n)
    n = re.sub(r"^[a-z]\d+_(?=[a-z]+z?\d*$)", "", n)
    return re.sub(r"_bb$", "", n)


def classify_channel(name, cohort=""):
    """
    Input: one channel name as written in the .vhdr, cohort (from cohort_of) for cohort overrides
    Returns: dictionary: type, subtype, group, contact, side, muscle, source (why we think so)
      lead: group lead1 (DBS/LFP 1-8) or lead2 (DBS 9-16, lab convention used by crosscheck.py), lead_left / lead_right for DBS1_L / DBS1_R; contact 1-8
      ecog: group ecog_anterior / ecog_posterior when the name says A/P, else ecog; contact = number in the name
      emg: muscle (full name, blank for emg_1), side (contra/ipsi/left/right, blank if not written)
      accel / joystick: contact = axis
    """
    n = normalize(name)
    out = {"type":"unknown", "subtype":"", "group":"", "contact":"", "side":"", "muscle":"", "source":"brainvision header name"}
    for pattern, kind, evidence in COHORT_OVERRIDES.get(cohort.upper(), []):
        if re.fullmatch(pattern, n):
            out.update(type=kind[0], subtype=kind[1], group=kind[1], source=f"cohort override: {evidence}")
            return out
    for pattern, kind in PATTERNS:
        m = re.fullmatch(pattern, n)
        if m is None:
            continue
        g = {k:v for k, v in m.groupdict().items() if v}
        out["type"], out["subtype"], out["group"] = kind[0], kind[1], kind[0]
        if kind == ("lfp", "lead"):
            k = int(g.get("contact", 0))
            if "hemi" in g:
                out["group"], out["contact"] = {"l":"lead_left", "r":"lead_right"}[g["hemi"]], str(k)
            elif k == 0: #"LFP" alone
                out["group"] = "lead"
            else:
                out["group"], out["contact"] = ("lead1", str(k)) if k <= 8 else ("lead2", str(k - 8))
            if k > 8 and "hemi" not in g:
                out["source"] = "brainvision header name; DBS9-16 = lead 2 contacts 1-8 by lab convention (crosscheck.lead_groups)"
        elif kind == ("lfp", "ecog"):
            out["group"] = {"a":"ecog_anterior", "p":"ecog_posterior", "":"ecog"}[g.get("a", g.get("b", ""))]
            out["contact"] = g["contact"]
        elif kind[0] == "emg":
            out["muscle"] = MUSCLES[g["muscle"]] + (f" ({g['part']})" if "part" in g else "")
            out["side"] = SIDES.get(g.get("side", ""), "")
        elif kind[0] == "mer" and kind[1] in ("micro", "macro"): #raw or spike band when the name says
            band = "_spike" if re.search(r"sp(i)?k", n) else "_raw" if "raw" in n else "_lfp" if "lfp" in n else ""
            out["subtype"] = kind[1] + band
        elif kind[0] in ("accel", "joystick"):
            out["contact"] = g["axis"]
        return out
    out["subtype"] = next((note for pattern, note in UNKNOWN_NOTES.items() if re.fullmatch(pattern, n)), "")
    return out


def channel_types(vhdr, cohort=None):
    """
    Input: a brainvision .vhdr, cohort (default: from the path)
    Returns: list of dictionaries, one per channel in file order: number (1-based), name, units, stream, plus everything from classify_channel.
             Channels whose units aren't volts (e.g. BIP, C, %) go to their own stream "<stream>_<units>" so nothing is labeled V that isn't.
    """
    cohort = cohort_of(vhdr) if cohort is None else cohort
    out = []
    channels = read_vhdr_channels(vhdr)
    names = [n.strip().lower() for n, _ in channels]
    for i, (name, units) in enumerate(channels):
        c = classify_channel(name, cohort)
        if names.count(name.strip().lower()) > 1: #BI002 12 Month Baseline: LFP1 on Ch20 and Ch25, which one is the contact isn't known
            c["source"] += f"; WARNING the header gives this name to {names.count(name.strip().lower())} channels"
        stream = STREAM_NAMES.get((c["type"], c["subtype"]), STREAM_NAMES.get(c["type"], "unknown"))
        if units not in VOLT_UNITS:
            stream = f"{stream}_{re.sub(r'[^A-Za-z0-9]+', '', units) or 'units'}"
        out.append({"number":i + 1, "name":name, "units":units, "stream":stream, **c})
    return out

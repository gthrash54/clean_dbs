import re

#one profile per cohort, keyed by the study ID prefix (lowercase). Only what the Box folder and file names show is written here;
#"project" is a best guess from those names until someone checks it against the protocol, so treat it as a hint, not a fact.
#layout: "stage" = <id>/stageN surgery/<recording> (uh, UG; output path <id>/stageN/<block>.h5, never change it)
#        "session" = recordings anywhere below <id>, the subfolder path is the session (output path <id>/<session>/<stem>.h5)
#modalities: rigs seen; project: what the study seems to be; folders: where recordings sit below the ID;
#names_seen: the evidence, recurring folder/file names (digits as #) that condition_of turns into condition labels; notes: quirks
COHORTS = {
    "uh": {"layout":"stage", "modalities":"tdt+bv", "project":"intraop DBS, stage 1 (TDT: ecog, lead, MER, EMG, audio, video) and stage 2 (BrainVision lead)",
           "folders":"<id>/stageN surgery/<block folder or .vhdr>", "names_seen":"dbs#, contact#, stn_dbs#, gpi_dbs#, baselinerest, reversal",
           "notes":"4 digit IDs (uh3021)"},
    "ug": {"layout":"stage", "modalities":"bv", "project":"intraop recordings, BrainVision in both stages: LFP baseline/pre/lesion, MER, macrostim, ERP assay, SD curve",
           "folders":"<id>/stageN surgery/<.vhdr>", "names_seen":"lfp_baseline, lfp_pre, lfp_lesion, mer, mer_stim, macrostim, erpassay, sdcurve, buzz",
           "notes":"4 digit IDs (UG3005); 'lesion' recordings suggest a lesioning case"},
    "dlpfc": {"layout":"session", "modalities":"bv", "project":"DLPFC stimulation with scalp EEG: rest (eyes open/closed), iTBS, paired pulses, bursts, cortical stim and AX-CPT, pre and post DBS",
              "folders":"files at <id>/, some in task data/axcpt/ieeg/, task data/motif protocol/, intraop/, preop/, eeg/",
              "names_seen":"baseline, rest, rest_eyesopen/closed, exp#_predbs/postdbs, itbs, pp, burstpp, corticalstim, #ma_probing, mepthreshold, practice, task#, simon",
              "notes":"~150 IDs, 3 or 4 digits, some lowercase folder names; BPP and AXCPT-### folders can sit inside an ID"},
    "bi": {"layout":"session", "modalities":"bv", "project":"postop follow-up stim and LFP: per-contact recordings, monopolar review, ARIE best/worst, home settings, capsular, ramp, ERPs",
           "folders":"<id>/LFP/, Stim/, ERPs/, Monopolar/, LEAP/, '# month(s)/' follow-up visits with their own LFP/, ramp/, home/, capsular/",
           "names_seen":"lfp#, stim#, contact#, arie_best/worst_#hz, home_settings_#hz, (pseudo)monopolar, capsular, ramp, teststim, rev",
           "notes":"follow-up visit comes from the '# month' folder (state field)"},
    "pp": {"layout":"session", "modalities":"bv", "project":"paired pulse stimulation, with ramp, sensory and resting LFP",
           "folders":"files at <id>/, LFP###/ and CHR###/ (other study IDs nested inside), OR_raw/, ramp/, SEP#/, NeuroOmega/",
           "names_seen":"pp#, ramp, sensory, baseline, lfp, stim#, #-# (contact pair)", "notes":"nested IDs are reported as their own patient with parent_id = the PP ID"},
    "mc": {"layout":"session", "modalities":"bv", "project":"bipolar contact-pair recordings, awake and asleep, with MER tests",
           "folders":"files at <id>/, Awake/, Asleep/", "names_seen":"#-# (contact pair), #-#_asleep, #-#_medtronic, mer_test, baseline",
           "notes":"awake/asleep comes from the folder or the file name (state field)"},
    "bv": {"layout":"session", "modalities":"bv", "project":"bipolar stim configurations (1+2-, 1-2+, ...) with EEG", "folders":"files at <id>/, HD_EEG/",
           "names_seen":"#+#-, #-#+ (stim polarity configuration)", "notes":""},
    "orbv": {"layout":"session", "modalities":"bv", "project":"bipolar stim configurations under anesthesia, as BV", "folders":"files at <id>/, ORBV#_Anesth/",
             "names_seen":"#+#-, #-#+, #hz", "notes":"Anesth folders are asleep"},
    "bpp": {"layout":"session", "modalities":"bv", "project":"bursts and paired pulses", "folders":"files at <id>/; some files sit straight in the number folder",
            "names_seen":"bpp#, burstsandpp", "notes":"also nested inside DLPFC IDs"},
    "bb": {"layout":"session", "modalities":"bv", "project":"LFP and MER with stim and bursts", "folders":"files at <id>/ (this is the BrainVision part only)",
           "names_seen":"lfp#, mer#, stim, buzz, rest, bursts, predbs/postdbs", "notes":""},
    "chr": {"layout":"session", "modalities":"bv", "project":"unclear from names, file names look like dates or contact sets", "folders":"files at <id>/; also nested in LFP and PP IDs",
            "names_seen":"#-#-#", "notes":""},
    "pt": {"layout":"session", "modalities":"bv", "project":"sensory testing, multipolar and directional stim, paired pulses", "folders":"files at <id>/; DIR###/, LFP###/, SEP##/ inside",
           "names_seen":"sensory, multipolar, lfp, pp, testing", "notes":"SEP## folders inside are nested IDs"},
    "sep": {"layout":"session", "modalities":"bv", "project":"sensory evoked potentials", "folders":"files at <id>/", "names_seen":"baseline", "notes":"2 digit IDs, some written SEP_##"},
    "pns": {"layout":"session", "modalities":"bv", "project":"stimulation with LFP (peripheral nerve stim going by the name)", "folders":"files at <id>/", "names_seen":"stim, lfp", "notes":""},
    "tm": {"layout":"session", "modalities":"bv", "project":"stim, paired pulses and bursts with LFP", "folders":"files at <id>/", "names_seen":"stim, lfp, pp, burst", "notes":""},
    "tmscs": {"layout":"session", "modalities":"bv", "project":"stim at set mA with rest recordings", "folders":"files at <id>/", "names_seen":"stim, #ma, rest", "notes":""},
    "lfp": {"layout":"session", "modalities":"bv", "project":"resting LFP and monopolar recordings", "folders":"files at <id>/; CHR### nested inside",
            "names_seen":"lfp, baseline, mono", "notes":"2 or 3 digit IDs; LFP### folders inside PP are this cohort, not sessions"},
    "joy": {"layout":"session", "modalities":"bv", "project":"joystick task with LFP and stim", "folders":"files at <id>/", "names_seen":"task (joystick), lfp, stim, pp, buzz", "notes":""},
    "ecog": {"layout":"session", "modalities":"bv", "project":"ECoG strip recordings", "folders":"files at <id>/", "names_seen":"ecog, #-#, traj, teststim", "notes":"2 digit IDs"},
    "burst": {"layout":"session", "modalities":"bv", "project":"burst stimulation", "folders":"files at <id>/", "names_seen":"", "notes":"also a BurstPairPulses folder with no ID"},
    "dys": {"layout":"session", "modalities":"bv", "project":"LFP (dystonia going by the name)", "folders":"files at <id>/", "names_seen":"lfp", "notes":""},
    "sensory": {"layout":"session", "modalities":"bv", "project":"sensory stim with LFP", "folders":"files at <id>/", "names_seen":"stim, lfp, initial", "notes":""},
    "pr": {"layout":"session", "modalities":"bv", "project":"polarity reversal (ARIE)", "folders":"files at <id>/; also a PR folder with no ID", "names_seen":"reversal, polarity", "notes":""},
    "motif": {"layout":"session", "modalities":"bv", "project":"DLPFC MOTIF protocol", "folders":"a MOTIF folder with no ID; DLPFC IDs also hold task data/motif protocol/",
              "names_seen":"#ma", "notes":""},
    "sdc": {"layout":"session", "modalities":"bv", "project":"strength-duration curve", "folders":"files at <id>/", "names_seen":"sdcurve", "notes":""},
    "or": {"layout":"session", "modalities":"mostly not BrainVision; a few IDs hold .vhdr", "project":"OR recordings; the BrainVision ones are bipolar stim configurations",
           "folders":"files at <id>/", "names_seen":"#+#-, #-#+", "notes":"54 OR IDs on Box, only the BrainVision part is covered here"},
}

#longest prefix first so tmscs isn't read as tm, orbv as or, bpp as bp
_prefixes = "|".join(sorted(COHORTS, key=len, reverse=True))
#a study ID is the prefix, an optional space/_/-, and 2-5 digits, nothing else ("PP1" or "ORBV1_Anesth" are sessions, not IDs)
ID_RE = re.compile(rf"^({_prefixes})[ _-]?(\d{{2,5}})$", re.IGNORECASE)
#an ID somewhere inside a file name, for recordings that sit straight in a number folder (e.g. "BPP012.vhdr", "dlpfc045_rest.vhdr")
ID_IN_NAME_RE = re.compile(rf"(?<![a-z])({_prefixes})[ _-]?(\d{{2,5}})(?!\d)", re.IGNORECASE)
#Box puts every person in a de-identified number folder (digits and dots, sometimes a suffix)
NUMBER_RE = re.compile(r"^\d+[.\d]*")


def cohort_of(study_id):
    """
    Input:
    study_id: a study ID folder name like "uh3021", "DLPFC045", "SEP_03"
    Returns: the cohort prefix in lowercase ("uh", "dlpfc", "sep"), or None if it isn't a study ID
    """
    match = ID_RE.match(str(study_id).strip())
    return match.group(1).lower() if match else None


#what was going on in a recording, read from its session folders and file name. first rule that matches wins, so specific before general.
#digits are kept as is, the patterns allow for them
CONDITIONS = [
    (r"imped[ae]nce", "impedance_check"), #impedance tests sit in task folders (DLPFC Task data/AXCPT/ieeg/impedencecheck), they aren't the task
    (r"axcpt", "task_axcpt"),
    (r"simon", "task_simon"),
    (r"motif", "task_motif"),
    (r"prac(tice)?", "task_practice"),
    (r"joystick", "task_joystick"),
    (r"(?<![a-z])(?<!no[ _])task", "task"), #"itbs no task" is not a task, same as "itbs_notask"
    (r"itbs", "stim_itbs"),
    (r"burst", "stim_burst"),
    (r"pairedpulse|(?<![a-z])b?pp(?![a-z])", "stim_paired_pulse"),
    (r"ramp", "stim_ramp"),
    (r"sdcurve|sdc(?![a-z])", "stim_sd_curve"),
    (r"erp", "stim_erp"),
    (r"mep ?thres?hold", "stim_mep_threshold"), #"MEPthrehold" is written on Box too
    (r"probing", "stim_probing"),
    (r"polarity|reversal", "stim_polarity_reversal"),
    (r"\d[+-]\d?[+-]", "stim_bipolar_config"),
    #programming-style recordings: named settings, contact configurations, or a frequency/amplitude in the name (BI "1_c_130hz_3ma")
    (r"arie|home_?settings|home_?pseudo|capsular|monopolar|pseudo|psuedo|multipolar|mono|ther(apeutic)?(?![a-z])|window|\d ?hz|\d ?ma(?![a-z])", "stim_setting"),
    (r"sensory|(?<![a-z])sep", "sensory"),
    (r"lesion", "lesion"),
    (r"buzz", "buzz"),
    (r"mer(?![a-z])|mer_", "mer"),
    (r"nostim", "rest"),
    (r"co[mn]tact", "contact"), #"comtact1" (BI032) is a typo for contact1
    (r"(?<!pre)stim|(?<![a-z])dbs\d", "stim"), #"restingstateprestimulation" is rest before stim, not stim
    (r"eyes ?open", "rest_eyes_open"),
    (r"eyes ?closed", "rest_eyes_closed"),
    (r"rest|baseline", "rest"),
    (r"lfp|ecog", "lfp"),
    (r"(?<![\d])\d+-\d+(?![\d-])", "contact_pair"),
]
#state of the patient or visit, separate from the condition, joined with ";" when several apply
STATES = [
    (r"asleep|anesth|precedex", "asleep"),
    (r"awake", "awake"),
    (r"pre_?dbs", "pre_dbs"),
    (r"post_?dbs", "post_dbs"),
    (r"pre_?op(?!en)", "preop"), #not "impedencetest_preopening" (before opening, intraop)
    (r"off_?meds", "off_meds"),
    (r"(?<![a-z])on_?meds", "on_meds"),
    (r"intraop", "intraop"),
    (r"(\d+) ?months?", "followup_{0}mo"),
]
#labels too vague to beat a more specific one from the folders
GENERIC = ("", "task", "stim", "contact", "rest", "lfp")
#the joystick cohort calls its task just "task"
COHORT_TASK = {"joy":"task_joystick"}


def condition_of(text, cohort=""):
    """
    Input:
    text: session folders and file name of one recording, e.g. "task data/axcpt/ieeg/task1" or "6 month/lfp/contact1_2"
    cohort: cohort prefix if known, for cohort-specific names
    Returns: (condition, state), "" when the name doesn't say. condition is one label (task_axcpt, stim_itbs, rest, ...),
        state is ";"-joined (asleep, pre_dbs, followup_6mo, ...)
    """
    text = str(text).lower().replace("\\", "/")
    folders, _, name = text.rpartition("/")
    #the file name usually says the most; a folder only wins when the name is blank or generic (task1 inside task data/axcpt/,
    #Bicep inside Initial Peripheral Stim/)
    condition, from_folders = _match(name), _match(folders)
    if (condition == "" and from_folders) or (condition in GENERIC and from_folders not in GENERIC):
        condition = from_folders
    if condition == "task" and cohort in COHORT_TASK:
        condition = COHORT_TASK[cohort]
    states = []
    for pattern, label in STATES:
        match = re.search(pattern, text)
        if match and label.format(*match.groups()) not in states:
            states.append(label.format(*match.groups()))
    return condition, ";".join(states)


def _match(text):
    return next((label for pattern, label in CONDITIONS if re.search(pattern, text)), "")

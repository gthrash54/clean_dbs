#version of the cleaned output, saved as the clean_version attribute at the top of every .h5
#bump it by 1 whenever a change to the code changes what goes in a cleaned file (new or renamed fields, different values),
#so clean.py re-cleans older files instead of skipping them; changes that leave the .h5 the same don't need a bump
#1: first numbered version (truncated .tev check, Chicago start time, stream-based duration, case_status, saturation at the rail)
#2: brainvision as typed streams (lead, ecog, emg, ...) with per-channel labels, read in windows; marker meanings, stim from markers.py
#   and /sidecars; stream rail attribute (brainvision actiCHamp too); brainvision start as the vmrk clock; cohort in info for every file;
#   CnA_ typed stim_params, PerB = voice B pulses, uneven tdt epocs trimmed (info epoc_warning); lead<N>_not_implanted and stage 2 labels
#3: lead contact coordinates per channel (channel_x/y/z/position, coord_space) from _cases/contacts.csv; analog_gain on tdt mer streams
CLEAN_VERSION = 3

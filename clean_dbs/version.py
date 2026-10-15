#version of the cleaned output, saved as the clean_version attribute at the top of every .h5
#bump it by 1 whenever a change to the code changes what goes in a cleaned file (new or renamed fields, different values),
#so clean.py re-cleans older files instead of skipping them; changes that leave the .h5 the same don't need a bump
#1: first numbered version (truncated .tev check, Chicago start time, stream-based duration, case_status, saturation at the rail)
#2: brainvision as typed streams (lead, ecog, emg, ...) with per-channel labels, read in windows; marker meanings, stim from markers.py
#   and /sidecars; stream rail attribute (brainvision actiCHamp too); brainvision start as the vmrk clock; cohort in info for every file;
#   CnA_ typed stim_params, PerB = voice B pulses, uneven tdt epocs trimmed (info epoc_warning); lead<N>_not_implanted and stage 2 labels
#3: lead contact coordinates per channel (channel_x/y/z/position, coord_space) from _cases/contacts.csv; analog_gain on tdt mer streams
#4: other cohorts, from a trial run on real recordings: condition/state rules (folder label when the name says nothing, the study ID in
#   a name isn't a condition, impedance_check, off/on_meds, typos), same ID nested in itself keeps its session; MC_Stimulus programs found
#   by their header, stricter sidecar matching, per-code matching, ambiguous TTL tables leave condition blank, TTL lists written as rows,
#   binary MC_Stimulus files listed; headers that repeat a channel name can be read
#5: inching blocks named without the underscore (uh3002 increment1depth3) get their depth's contact coordinates, not the final ones
#6: an implanted lead never plugged into the amplifier (cases.csv note "<target> lead never connected", uh3016 GPi) is labeled
#   lead<N>_not_connected with no target
#7: per-channel channel_dnu (1 = the lab marked it do not use: bad contacts, contacts still in the cannula) and stream attribute dnu_list,
#   from _cases/dnu_channels.csv (scripts/make_dnu_table.py builds it from the lab's "UH3 or notes - data to not use.xlsx")
CLEAN_VERSION = 7

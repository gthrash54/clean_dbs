function out = clean_dbs_load(h5file, varargin)
%CLEAN_DBS_LOAD Loads a cleaned clean_dbs .h5 file into one MATLAB struct (no Python needed).
%   out = clean_dbs_load(h5file)
%   out = clean_dbs_load(h5file, 'types', {'lfp'}, 't1', 10, 't2', 20)
%
%   Input:
%   h5file: a cleaned .h5 file made by scripts/clean.py (e.g. <out>/uh3021/stage1/clinical_dbs1.h5)
%   optional name/value pairs:
%     'types': cell of stream types to load, e.g. {'lfp'} or {'mer','emg'} (default: every stream)
%              types are lfp, mer, emg, audio, stim_monitor, unknown; epocs, scalars and flags always load
%     't1', 't2': time window in seconds from the start of the recording (default 0 and 0 = whole recording)
%              only the samples inside the window are read from disk, so this is fast even on 1+ GB files
%     'time': true adds a .time vector (seconds) to each stream (default false, it costs as much memory as a channel)
%   Returns: out, a struct laid out like the file
%     out.streams.<name>: data (channels x samples, single, volts), fs, start_time (time of the first sample read),
%                         channels, device, type, subtype, units, name (the name in the file)
%     out.epocs.<name>:   onset, offset, data (TDT only), type, name; only events with t1 <= onset < t2
%     out.scalars.<name>: ts, data (parameters x events, e.g. eS1p stim settings), type, name; only t1 <= ts < t2
%     out.info:           block, start, duration, subject, video, ... (every attribute the reader saved)
%     out.flags:          saturation_<stream> (channel, start, end), new_segments (time, missing, is_gap), dnu
%                         flags always cover the whole recording, even with a time window
%     out.provenance:     source, clean_dbs_commit, created, package versions
%   Names that aren't valid MATLAB field names are fixed with matlab.lang.makeValidName; .name keeps the original.
%   Numbers come back as double, except stream data which stays single (double would double the memory).
%
%   Example:
%     r = clean_dbs_load('clinical_dbs1.h5', 'types', {'lfp'}, 't1', 10, 't2', 20);
%     x = r.streams.ecos.data;                                   % 28 x samples, volts
%     t = r.streams.ecos.start_time + (0:size(x,2)-1) / r.streams.ecos.fs;
%   Contact Garrett Thrash for any issues garrettthrash54@gmail.com

p = inputParser;
p.addRequired('h5file', @(x) ischar(x) || isstring(x));
p.addParameter('types', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('t1', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('t2', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('time', false, @(x) islogical(x) || isnumeric(x));
p.parse(h5file, varargin{:});
h5file = char(p.Results.h5file);
types = cellstr(p.Results.types); %so 'lfp' and {'lfp'} both work
t1 = double(p.Results.t1);
t2 = double(p.Results.t2);
if ~isfile(h5file)
    error('clean_dbs_load:notFound', 'File not found: %s, check path and spelling', h5file);
end
if t1 < 0 || (t2 ~= 0 && t2 <= t1)
    error('clean_dbs_load:badWindow', 'Bad time window t1=%g, t2=%g: need t1 >= 0 and t2 = 0 (to the end) or t2 > t1', t1, t2);
end
if t2 == 0
    t2 = inf; %0 means read to the end, same as read_tdt / read_bv
end

meta = h5info(h5file); %only the layout and attributes, no data is read yet
out = struct('streams', struct(), 'epocs', struct(), 'scalars', struct(), 'info', struct(), 'flags', struct(), 'provenance', attrs(meta));

%streams: read only the samples inside the window
for g = child_groups(meta, '/streams')
    s = attrs(g);
    if ~isempty(types) && ~any(strcmp(s.type, types))
        continue
    end
    s.name = leaf(g.Name);
    s.channels = read_ds(h5file, g, 'channels');
    for d = {g.Datasets.Name} %what each channel is (channel_names, channel_group, channel_target, ...), as text
        if startsWith(d{1}, 'channel_')
            s.(d{1}) = cellstr(string(h5read(h5file, [g.Name '/' d{1}])));
        end
    end
    ds = find_ds(g, 'data');
    n_total = ds.Dataspace.Size(1); %h5py writes channels x samples, MATLAB sees it flipped: samples x channels
    first = max(0, ceil((t1 - s.start_time) * s.fs - 1e-9)); %first sample at or after t1 (0-based)
    last = min(n_total, ceil((t2 - s.start_time) * s.fs - 1e-9)); %one past the last sample before t2
    n = last - first;
    if n > 0
        s.data = h5read(h5file, [g.Name '/data'], [first + 1, 1], [n, numel(s.channels)])'; %transpose back to channels x samples
    else
        s.data = zeros(numel(s.channels), 0, 'single');
    end
    s.start_time = s.start_time + first / s.fs; %time of the first sample we actually read
    if p.Results.time
        s.time = s.start_time + (0:size(s.data, 2) - 1) / s.fs;
    end
    out.streams.(valid(s.name)) = s;
end

%epocs and scalars are small, so read them whole and keep the events inside the window
for g = child_groups(meta, '/epocs')
    e = attrs(g);
    e.name = leaf(g.Name);
    for d = {g.Datasets.Name}
        e.(d{1}) = read_ds(h5file, g, d{1});
    end
    keep = e.onset >= t1 & e.onset < t2;
    for d = {'onset', 'offset', 'data'}
        if isfield(e, d{1})
            e.(d{1}) = e.(d{1})(keep);
        end
    end
    out.epocs.(valid(e.name)) = e;
end
for g = child_groups(meta, '/scalars')
    sc = attrs(g);
    sc.name = leaf(g.Name);
    sc.ts = read_ds(h5file, g, 'ts');
    sc.data = read_ds(h5file, g, 'data'); %parameters x events
    keep = sc.ts >= t1 & sc.ts < t2;
    sc.ts = sc.ts(keep);
    sc.data = sc.data(:, keep);
    out.scalars.(valid(sc.name)) = sc;
end

%info is all attributes; flags are attributes (dnu) plus one group per flag
info = child_groups(meta, '');
info = info(strcmp({info.Name}, '/info'));
if ~isempty(info)
    out.info = attrs(info);
end
flags = child_groups(meta, '');
flags = flags(strcmp({flags.Name}, '/flags'));
if ~isempty(flags)
    out.flags = attrs(flags);
    for g = child_groups(meta, '/flags')
        f = struct();
        for d = {g.Datasets.Name}
            f.(valid(d{1})) = read_ds(h5file, g, d{1});
        end
        out.flags.(valid(leaf(g.Name))) = f;
    end
end
end


function groups = child_groups(meta, path)
%Input: h5info of the file, a group path like '/streams' ('' = top level). Returns: its child groups as a row (empty if missing)
groups = meta.Groups;
if ~isempty(path)
    groups = groups(strcmp({groups.Name}, path));
    if isempty(groups)
        return
    end
    groups = groups.Groups;
end
groups = reshape(groups, 1, []); %row, so "for g = groups" loops one group at a time
end


function s = attrs(g)
%Input: one group from h5info. Returns: struct of its attributes, text as char, numbers as double
s = struct();
for a = reshape(g.Attributes, 1, [])
    v = a.Value;
    if iscell(v) && isscalar(v) %some MATLAB versions wrap h5py text in a cell
        v = v{1};
    end
    if isstring(v)
        v = char(v);
    elseif isnumeric(v)
        v = double(v); %flags like dnu are saved as int8, double is easier to work with
    end
    s.(valid(a.Name)) = v;
end
end


function x = read_ds(h5file, g, name)
%Input: file, group from h5info, dataset name. Returns: the dataset as double, 2-D ones flipped back to the Python layout
ds = find_ds(g, name);
if isempty(ds.Dataspace.Size) || any(ds.Dataspace.Size == 0)
    x = zeros(0, 1); %h5read errors on empty datasets (e.g. no saturation found)
    return
end
x = double(h5read(h5file, [g.Name '/' name]));
if numel(ds.Dataspace.Size) == 2
    x = x'; %h5py rows x columns comes in as columns x rows (checked on the file's shape, so 1 x events stays 1 x events)
end
end


function ds = find_ds(g, name)
ds = g.Datasets(strcmp({g.Datasets.Name}, name));
end


function name = leaf(path)
%'/streams/ecos' -> 'ecos'
parts = strsplit(path, '/');
name = parts{end};
end


function name = valid(name)
%names like 'mic_' are fine, but anything starting with a number or holding odd characters isn't a legal field name
name = matlab.lang.makeValidName(name);
end

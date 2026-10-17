function ID = bebo_ids(who)
% BEBO_IDS  Participant IDs of the analysed sample (93: 37 YA, 56 OA).
%   who: 'ya' | 'oa' | 'all' (YA first, then OA -- the order used by every
%   group-level step). Must match subj_list_all93.txt.
% Excluded: 3111, 3213 (YA; overall proportion correct .20/.21, below
% chance and below every included participant); 4210 (OA; excluded in the
% original sample).
YA = [3109,3110,3112,3113,3114,3115,3117,3118,3119,3120,3121,3122,3123,3124,3125,3126,3127,3128,3129,3130,3131,3132, ...
      3205,3206,3207,3208,3209,3211,3212,3214,3215,3216,3217,3218,3219,3220,3222];
OA = [4109,4110,4111,4112,4113,4114,4115,4116,4118,4119,4120,4121,4122,4123,4124,4125,4126,4127,4128,4129,4131,4132,4133, ...
      4134,4135,4136,4137,4138,4207,4208,4209,4211,4212,4213,4214,4215,4216,4217,4218,4219,4220,4221,4222,4223,4224,4225, ...
      4227,4228,4229,4230,4231,4232,4233,4236,4238,4239];
if ~isempty(getenv('BEBO_IDS_FILE'))   % test runs only: line 1 = YA IDs, line 2 = OA IDs
    L = strsplit(strtrim(fileread(getenv('BEBO_IDS_FILE'))), newline);
    YA = str2num(L{1}); OA = str2num(L{2}); %#ok<ST2NM>
end
switch lower(who)
    case 'ya',  ID = YA;
    case 'oa',  ID = OA;
    case 'all', ID = [YA, OA];
    otherwise,  error('bebo_ids:who', 'Unknown group "%s" (ya, oa or all)', who)
end
end

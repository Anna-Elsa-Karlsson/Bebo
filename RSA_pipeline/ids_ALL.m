function [ID] = ids_ALL(who,varargin)
% BEBO IDs
% Input:
% - who (STRING): ya, ch, all 
% Output:
% - ID (DOUBLE)

if nargin == 1
    % 2026-10-19 (user): ALL 96 with usable EEG -- 3111, 3213, 4210 back in.
    % They were removed on 2026-09-29 only because no pre-made TFR files
    % existed for them; the rerun decomposes everything from PreprocData.
    YA = [3109,3110,3111,3112,3113,3114,3115,3117,3118,3119,3120,3121,3122,3123,3124,3125,3126,3127,3128,3129,3130,3131,3132, ...
         3205,3206,3207,3208,3209,3211,3212,3213,3214,3215,3216,3217,3218,3219,3220,3222];

%     CH = [2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108, 2109, 2110, 2111, 2112, 2113, 2114, 2115, 2116, 2201, 2202, 2203, 2204, 2205, 2206];
%
    OA = [4109,4110,4111,4112,4113,4114,4115,4116,4118,4119,4120,4121,4122,4123,4124,4125,4126,4127,4128,4129,4131,4132,4133, ...
        4134,4135,4136,4137,4138,4207,4208,4209,4210,4211,4212,4213,4214,4215,4216,4217,4218,4219,4220,4221,4222,4223,4224,4225, ...
        4227,4228,4229,4230,4231,4232,4233,4236,4238,4239];
    
elseif nargin > 1 && strcmp(varargin{1},'tutorial')
    YA = [3113,3115,3119,3122,3132,3208,3213,3214,3218,3219];
     
    CH = [2102,2104,2110,2112,2113,2114,2201,2203,2205,2206];
    % CH = {'02CSS25', '02DSW27', '03JBJ25',  '04CJC10',  '06KRB20', '07KNJ16',  '09KRE29',  '10SSC13', '11SWC30', '12SLZ15'}; % old names
end

% which IDs
if strcmp(who,'all') % all IDs for final analyses        
    ID = [YA, OA]; %CH, 
elseif strcmp(who,'ch') || strcmp(who,'CH') 
    ID = CH;
elseif strcmp(who,'ya') || strcmp(who,'YA')
    ID = YA;
elseif strcmp(who,'oa') || strcmp(who,'OA')
    ID = OA;
end
    
    
end

    

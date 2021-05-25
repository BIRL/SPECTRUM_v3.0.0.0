%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SPECTRUM: A MATLAB Toolbox for Top-down Proteomics     %
%                           Version 3.0.0.0                        %
%        Copyright (c) Biomedical Informatics Research Laboratory, %
%          Lahore University of Management Sciences Lahore (LUMS), %
%                           Pakistan.                              %
%                (http://biolabs.lums.edu.pk/BIRL)                 %
%                    (safee.ullah@gmail.com)                       %
%                 Last Modified on: 25-May-2021                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fix_lines(fname, fname2)
%FIX_LINES  Improves the line style of eps files generated by print
%
% Examples:
%   fix_lines fname
%   fix_lines fname fname2
%
% This function improves the style of lines in eps files generated by
% MATLAB's print function, making them more similar to those seen on
% screen. Grid lines are also changed from a dashed style to a dotted
% style, for greater differentiation from dashed lines.
%
% The function also places embedded fonts after the postscript header, in
% versions of MATLAB which place the fonts first (R2006b and earlier), in
% order to allow programs such as Ghostscript to find the bounding box
% information.
%
% IN:
%   fname - Name or path of source eps file.
%   fname2 - Name or path of destination eps file. Default: same as fname.

% Copyright: (C) Oliver Woodford, 2008-2010

% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928)
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)

% Thank you to Sylvain Favrot for bringing the embedded font/bounding box
% interaction in older versions of MATLAB to my attention.
% Thank you to D Ko for bringing an error with eps files with tiff previews
% to my attention.
% Thank you to Laurence K for suggesting the check to see if the file was
% opened.

% Read in the file
fh = fopen(fname, 'r');
if fh == -1
    error('File %s not found.', fname);
end
try
    fstrm = fread(fh, '*char')';
catch ex
    fclose(fh);
    rethrow(ex);
end
fclose(fh);

% Move any embedded fonts after the postscript header
if strcmp(fstrm(1:15), '%!PS-AdobeFont-')
    % Find the start and end of the header
    ind = regexp(fstrm, '[\n\r]%!PS-Adobe-');
    [ind2 ind2] = regexp(fstrm, '[\n\r]%%EndComments[\n\r]+');
    % Put the header first
    if ~isempty(ind) && ~isempty(ind2) && ind(1) < ind2(1)
        fstrm = fstrm([ind(1)+1:ind2(1) 1:ind(1) ind2(1)+1:end]);
    end
end

% Make sure all line width commands come before the line style definitions,
% so that dash lengths can be based on the correct widths
% Find all line style sections
ind = [regexp(fstrm, '[\n\r]SO[\n\r]'),... % This needs to be here even though it doesn't have dots/dashes!
    regexp(fstrm, '[\n\r]DO[\n\r]'),...
    regexp(fstrm, '[\n\r]DA[\n\r]'),...
    regexp(fstrm, '[\n\r]DD[\n\r]')];
ind = sort(ind);
% Find line width commands
[ind2 ind3] = regexp(fstrm, '[\n\r]\d* w[\n\r]');
% Go through each line style section and swap with any line width commands
% near by
b = 1;
m = numel(ind);
n = numel(ind2);
for a = 1:m
    % Go forwards width commands until we pass the current line style
    while b <= n && ind2(b) < ind(a)
        b = b + 1;
    end
    if b > n
        % No more width commands
        break;
    end
    % Check we haven't gone past another line style (including SO!)
    if a < m && ind2(b) > ind(a+1)
        continue;
    end
    % Are the commands close enough to be confident we can swap them?
    if (ind2(b) - ind(a)) > 8
        continue;
    end
    % Move the line style command below the line width command
    fstrm(ind(a)+1:ind3(b)) = [fstrm(ind(a)+4:ind3(b)) fstrm(ind(a)+1:ind(a)+3)];
    b = b + 1;
end

% Find any grid line definitions and change to GR format
% Find the DO sections again as they may have moved
ind = int32(regexp(fstrm, '[\n\r]DO[\n\r]'));
if ~isempty(ind)
    % Find all occurrences of what are believed to be axes and grid lines
    ind2 = int32(regexp(fstrm, '[\n\r] *\d* *\d* *mt *\d* *\d* *L[\n\r]'));
    if ~isempty(ind2)
        % Now see which DO sections come just before axes and grid lines
        ind2 = repmat(ind2', [1 numel(ind)]) - repmat(ind, [numel(ind2) 1]);
        ind2 = any(ind2 > 0 & ind2 < 12); % 12 chars seems about right
        ind = ind(ind2);
        % Change any regions we believe to be grid lines to GR
        fstrm(ind+1) = 'G';
        fstrm(ind+2) = 'R';
    end
end

% Isolate line style definition section
first_sec = strfind(fstrm, '% line types:');
[second_sec remaining] = strtok(fstrm(first_sec+1:end), '/');
[remaining remaining] = strtok(remaining, '%');

% Define the new styles, including the new GR format
% Dot and dash lengths have two parts: a constant amount plus a line width
% variable amount. The constant amount comes after dpi2point, and the
% variable amount comes after currentlinewidth. If you want to change
% dot/dash lengths for a one particular line style only, edit the numbers
% in the /DO (dotted lines), /DA (dashed lines), /DD (dot dash lines) and
% /GR (grid lines) lines for the style you want to change.
new_style = {'/dom { dpi2point 1 currentlinewidth 0.08 mul add mul mul } bdef',... % Dot length macro based on line width
    '/dam { dpi2point 2 currentlinewidth 0.04 mul add mul mul } bdef',... % Dash length macro based on line width
    '/SO { [] 0 setdash 0 setlinecap } bdef',... % Solid lines
    '/DO { [1 dom 1.2 dom] 0 setdash 0 setlinecap } bdef',... % Dotted lines
    '/DA { [4 dam 1.5 dam] 0 setdash 0 setlinecap } bdef',... % Dashed lines
    '/DD { [1 dom 1.2 dom 4 dam 1.2 dom] 0 setdash 0 setlinecap } bdef',... % Dot dash lines
    '/GR { [0 dpi2point mul 4 dpi2point mul] 0 setdash 1 setlinecap } bdef'}; % Grid lines - dot spacing remains constant

if nargin < 2
    % Overwrite the input file
    fname2 = fname;
end

% Save the file with the section replaced
fh = fopen(fname2, 'w');
if fh == -1
    error('Unable to open %s for writing.', fname2);
end
try
    fwrite(fh, fstrm(1:first_sec), 'char*1');
    fwrite(fh, second_sec, 'char*1');
    fprintf(fh, '%s\r', new_style{:});
    fwrite(fh, remaining, 'char*1');
catch ex
    fclose(fh);
    rethrow(ex);
end
fclose(fh);
return
% Defines sl for separating folders
sl = filesep; % File folder seperator (sl) different on Mac and PC
AnalysesFolder = mfilename('fullpath');
TaskFolder = AnalysesFolder(1:max(strfind(AnalysesFolder,[sl 'Analyses'])));
MonkeyNames = {'Woody' 'Raul'};
SummarySessions = struct('Woody',[], 'Raul',[]);
for monkey = 1:length(MonkeyNames)
    
    monkeySessionCounter = 0;
    MonkeyName = MonkeyNames{monkey};
    
    % Get sessions
    ResultsFolders = dir(...
        [TaskFolder 'Results' sl MonkeyName sl '2*']);
    sessList = {ResultsFolders.name}';
    
    for sess = 1:length(sessList)
        
        SessionDate = sessList{sess};
        SessionFolder = [TaskFolder 'Results' sl MonkeyName sl SessionDate];
        files = dir(SessionFolder);
        dirFiles = {files.name};
        EyeCalFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'EyeCal')))));
    EyeEndFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'EyeEnd')))));
    XMazeFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'XMaze')))));
    FreeRoamFile = cell2mat(dirFiles(~cellfun(@isempty,(strfind(dirFiles,'Free')))));

    cd(SessionFolder)
    SummarySessions.(MonkeyName)
    if ~isempty(XMazeFile)
        load([SessionFolder sl XMazeFile]);
        length(fieldnames(XMazeStruct.Trials))
       XMazeStruct = Remove1000hzsampling(XMazeStruct,1);
        
        save(XMazeFile,'XMazeStruct')

    else display('no XMazeData');XMazeStruct = [];
    end
     if ~isempty(EyeEndFile)
        load([SessionFolder sl EyeEndFile]);
        EyeEndStruct = Remove1000hzsampling(EyeEndStruct,0);
                 save(EyeEndFile,'EyeEndStruct')

    else display('no EyeEnd data'); EyeEndStruct = [];
    end
    if ~isempty(EyeCalFile);
        load([SessionFolder sl EyeCalFile]);
       EyeCalStruct = Remove1000hzsampling(EyeCalStruct,0);
         save(EyeCalFile,'EyeCalStruct')
    else display('no EyeCal data'); EyeCalStruct = [];
    end
    if ~isempty(FreeRoamFile)
        load([SessionFolder sl FreeRoamFile]);
        FreeRoamStruct = Remove1000hzsampling(FreeRoamStruct,1);
        save(FreeRoamFile,'FreeRoamStruct')

    else display('no FreeRoam data'); FreeRoamStruct = [];
    end

    end% session loop
end% monkey loop
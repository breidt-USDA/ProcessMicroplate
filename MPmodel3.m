classdef MPmodel3 < handle
   % Matlab MPmodel3 model class:
   % Version: 3.240607
   % 
   % Author:
   % Fred Breidt, USDA/ARS
   % USDA/ARS Research Microbiologist
   % Food Science and Market Quality and Handling Research Unit
   % 322 Schaub Hall, Box 7624
   % NC State University, Raleigh, NC 27695-7624
   % Email: fred.breidt@usda.gov
   %
   % Program Description:
   % Process 96 well microtiter plate growth curves from optical density
   % readings generated using an automated microtiter plate reader.
   %
   % MPmodel3 class has functions for processing microtiter plate data to
   % determing microbial growth parameters. Parameters include growth rate 
   % (GR), lag time (LAG), initial or minimum OD (MIN), and maximum OD 
   % (MAX). The paramters are obtained using a sliding window method 
   % (Breidt, 1994). 
   % 
   % The calculated growth parameters (GR, LAG and MAX) are then used to
   % estimate a predicted Gompertz curve based on a published algorithm 
   % (Zewitering et al., 1992).   
   %
   % The optical density data used for calculation of the growth parameters
   % is modified by subtracting a background value from each OD reading.
   % The background value can be automatically obtained from selected
   % row(s)-column(s) values or manually input. 
   %
   % Data can be output from the program to .csv tables including a table
   % with growth parameters from selected row(s)-col(s), a stats table with
   % mean and standard deviations for replicate data, and a table with the
   % growth curve OD values used for determining the parameters.
   %
   % References:
   % 1. Breidt F, Romick TL, Fleming HP. 1994. A rapid method for the 
   % determination of bacterial growth kinetics. J Rapid Methods Autom 
   % Microbiol 3(1):59-68. 
   % https://doi.org/10.1111/j.1745-4581.1994.tb00308.x
   %   
   % 2. Zwietering MH, Jongenburger I, Rombouts FM, van't Tiet K. 1990. 
   % Modeling the bacterial growth curve. Appl Env Microbiol 
   % 56(6):1875-1881. https://doi.org/10.1128/aem.56.6.1875-1881.1990
   %
   % Public access functions:
   %     1) Constructor: MPmodel3(MP_XYY, dataname)
   %        Private: (internal data structures)
   %           MP_XYY      %Microplate data matrix 
   %           MP_XYN      %Microplate in 8x12xN data matrix format 
   %           TPs         %number of time points (including time zero)
   %        Public: (variables that can be accessed/modified)
   %           MachineRes  %resolution of spectrophotometer
   %           FirstN      %initial points (1..N) for determining minOD val
   %           blankData   %struct with blank data for subtracting from OD
   %              (see CrvData.blankdata, below for data struct members)
   %              NOTE: blankData.mean is used to subtract from OD values
   %                 for calculation of growth rate parameters
   %           dataname    %string with name of data object
   %        Returns:
   %           class object
   %     2) curvedata = getCurveData(rowvec, colvec, blankrowvec, 
   %                    blankcolvec, processbyN)
   %           rowvec = one or more row value(s)
   %           colvec = one or more col value(s)
   %           blankrowvec = one or more row values for blank val(s)
   %           blankcolvec = one or more col values for blank val(s)
   %           processbyN = size of sliding window for sets of data points
   %        Returns:
   %           CrvData.XYY = XYY;                  %assign XYY matrix
   %           CrvData.XY = XY;                    %assign XY
   %           CrvData.LnXY = temp.ObsLnXY;        %Ln of Nx2 data
   %           CrvData.PredLnXY = temp.PredLnXY;   %XY for Gompertz model
   %           CrvData.RowColTable:
   %              variables:
   %                 ROW','COL','GR','DBL','LAG','MIN','MAX'
   %              NOTE: MIN and LAG are zero if MIN OD value is at machine 
   %              resolution. Also Min and Max OD values are real numbers,
   %              not natural log, which is used to calculate growth rate.          
   %           CrvData.StatsTable, 
   %              variables: 
   %                 'GR','DBL','LAG','MIN','MAX'
   %              row names:
   %                 1. 'Mean', 2. 'Stdev'
   %           CrvData.blankdata,
   %              blankdata.XY = blank values (y = OD) in xy matrix
   %              blankdata.Nvals = number of OD values used
   %              blankdata.mean = mean of OD values
   %              blankdata.median = median of OD values
   %              blankdata.mode = mode of OD values
   %              blankdata.std = standard deviation of OD values
   %  Display functions: curvedata input is return struct from getCurveData
   %     3) showBlank() - histogram of all current blank values
   %     4) showMinOD(curvedata,cutoff) - bar graph, minOD for each row-col
   %              NOTE: cutoff is horizontal line at arbitrary minimum val
   %     5) showGR(curvedata) - bar graph of GR for each row-col 
   %     6) showLag(curvedata) - bar graph of lag values for each row-col
   %     7) showMaxOD(curvedata) - bar graph of maxOD for each row-col
   %  Write data to files:
   %     8) writeRowColTable(obj,rowvec,colvec,blankrowvec, ...
   %           blankcolvec, processbyN)
   %     9) writeCurveTable(obj,rowvec,colvec,blankrowvec, ...
   %           blankcolvec, processbyN)
   % Dependencies:
   %    LinReg.m function (linear regression algorithm  
   %----------------------------------------------------------------------
   properties (Access = private)    
      MP_XYY      %Microplate data matrix 
      MP_XYN      %Microplate in 8x12xN data matrix format 
      TPs         %number of time points (including time zero)
   end% private properties

   properties %public
      MachineRes  %resolution of spectrophotometer
      FirstN      %initial points (1..N) for determining min OD val
      blankData   %info about the blank value used for subtracting from OD
      dataname    %string with name of data object
      version     %version number 
   end% properties

   methods %public methods------------------------------------------------
      function obj = MPmodel3(MP_XYYtable, dataname)
      %constructor
         %take Nx97 Table and make MP_XYY matrix, MP_XYN 3D matrix, etc. 
         obj.MP_XYY = MP_XYYtable{:,:};      %save original data (8 x 97)
         obj.TPs = size(obj.MP_XYY,1);       %number of data points
         obj.MP_XYN = zeros(8,12,obj.TPs);   %create XYN 3D matrix
         obj.setXYN(0);                      %fill MP_XYN matrix
         obj.MachineRes = 0.001;             %machine resolution
         obj.FirstN = 3;                     %initial datapts for min OD
         obj.blankData = struct([]);         %blank values matrix Nx2
         obj.dataname = dataname;            %set data name
         obj.version = 3.240606;             %version (version.date-format)
      end %function

      function CrvData = getCurveData(obj, rowvec, colvec, blankrowvec, ...
            blankcolvec, processbyN)
      %get curves and GR data from one growth curve
         %set blank data struct (obj.blankData structure)
         blankdata = obj.setBlankData(blankrowvec,blankcolvec);
         rowcolmat = obj.getRowColMat(rowvec, colvec); %convert to matrix
         temp = obj.getGRdata(rowcolmat,processbyN); %get GR data struct
         XYY = obj.getSubXYY(rowcolmat);     %NxMultipleYs
         XY = obj.XYY2XY(XYY);               %Nx2 with reps (if present)
         CrvData.XYY = XYY;                  %assign XYY matrix
         CrvData.XY = XY;                    %assign XY
         CrvData.LnXY = temp.ObsLnXY;        %Ln of Nx2 data
         CrvData.PredLnXY = temp.PredLnXY;   %XY for Gompertz model
         %return tables with summarized data, convert Min and Max to real
         CrvData.RowColTable = obj.getResultsTable(rowcolmat,processbyN);
         %convert Ln OD values to real numbers in row col table 
         CrvData.RowColTable.MIN = exp(CrvData.RowColTable.MIN);
         %if MIN = log(MachineRes), show as zero, but correct for rounding
         %CrvData.RowColTable.MIN(CrvData.RowColTable.MIN<=(0.001+1e-12))=0;
         CrvData.RowColTable.MIN(CrvData.RowColTable.MIN<=0.001)=0;
         %convert Ln to real
         CrvData.RowColTable.MAX = exp(CrvData.RowColTable.MAX);
         %set statistics table
         CrvData.StatsTable = obj.getStatsTable(CrvData.RowColTable);
         CrvData.blankData = blankdata;      %assign blank value data
      end %function

      function writeRowColTable(obj,rowvec,colvec,blankrowvec, ...
         blankcolvec, processbyN)
      %write RowColTable data to a file or append existing one if in exists
         curvedata = obj.getCurveData(rowvec, colvec, blankrowvec, ...
            blankcolvec, processbyN);        %generate curve data 
         currentpath = pwd;                  %get current path
         datafile = append(obj.dataname,"_data",".csv"); %StatsTable name
         datapathname = fullfile(currentpath, datafile); %StatsTable path
         try
            writetable(curvedata.RowColTable,datapathname, ...
               'WriteMode','append');        %add data to file
         catch
            msg = "File(s) could not be written, " + ...
               "check pathname and/or " + ...
               "try closing .csv file";      %suggest problem fix
            error(msg);                      %post error message
         end %try-catch      end %function
      end %function

      function writeCurveTable(obj,rowvec,colvec, ...
         blankrowvec, blankcolvec, processbyN)
      %write curve table with XY data, LnXY and predicted LnXY
      %this function does not append data and will overwrite 
      % an existing file.
         curvedata = obj.getCurveData(rowvec, colvec, blankrowvec, ...
            blankcolvec, processbyN);        %generate curve data 
         currentpath = pwd;                  %get current path
         crvfile = append(obj.dataname,"_curve",".csv"); %curve table name
         crvpathname = fullfile(currentpath, crvfile); %curve table path
         obspred = [curvedata.XY, curvedata.LnXY(:,2), ...
            curvedata.PredLnXY(:,2)];        %make matrix for curve data
         CrvTable = array2table(obspred,'Variablenames',{'X','Y','LnY', ...
            'PredLnY'});                     %convert matrix to a table
         try
            writetable(CrvTable,crvpathname); %write data to file
         catch
            msg = "File(s) could not be written, " + ...
               "check pathname and/or " + ...
               "try closing .csv file";      %suggest problem fix
            error(msg);                      %post error message
         end %try-catch
      end %function

      function statpathname = writeStatsTable(obj,rowvec,colvec,blankrowvec, ...
         blankcolvec, processbyN,fileID)
      %write stat file or append to an existing one if it exists
         curvedata = obj.getCurveData(rowvec, colvec, blankrowvec, ...
            blankcolvec, processbyN);        %generate curve data 
         currentpath = pwd;                  %get current path
         statfile = append(obj.dataname,"_",fileID,".csv"); %StatsTable name
         statpathname = fullfile(currentpath, statfile); %StatsTable path
         stattable = curvedata.StatsTable;   %assign statistics table name
         try
            writetable(stattable,statpathname, ...
               'WriteRowNames',true,'WriteMode','append'); %add data 
         catch
            msg = "File(s) could not be written, " + ...
               "check pathname and/or " + ...
               "try closing .csv file";      %suggest problem fix
            error(msg);                      %post error message
         end %try-catch
      end %function

      function blankdata = setBlankData(obj,rowvec,colvec)
      %set blank value and return a struct with blank data stats
         rowcolmat = obj.getRowColMat(rowvec, colvec); %convert to matrix
         XYY = obj.getSubXYY(rowcolmat);
         XY = obj.XYY2XY(XYY);
         blankdata.XY = XY;
         blankdata.Nvals = size(XY,1);
         blankdata.mean = mean(XY(:,2));
         blankdata.median = median(XY(:,2));
         blankdata.mode = mode(XY(:,2));
         blankdata.std = std(XY(:,2));
         obj.blankData = blankdata; 
      end %function

      function Bh = showBlankData(obj)
      %plot histogram of blank values using blankData struct if it exists
        if ~isempty(fieldnames(obj.blankData)) %check for blank data
            set(gcf, 'NextPlot', 'new'); % next plot goes in new figure
           Bh = histogram(obj.blankData.XY(:,2)); %plot histogram 
        end
      end %function 

      function Mh = showMinOD(obj, curvedata, cutoffOD)
      %show graph of minimum values from RoCol table curve data 
         set(gcf, 'NextPlot', 'new');        %next plot goes in new figure
         RCtab = curvedata.RowColTable;      %get GR data all rows/cols
         RCtab = sortrows(RCtab,1);          %sort table by row vals    
         initvec = RCtab.MIN;                %extract minOD vec (as real)
         Xlabelvec = obj.getLabelVec([RCtab.ROW, RCtab.COL]); %Xaxis labels
         Mh = bar(Xlabelvec, initvec);       %bar graph of min with labels
         yline(cutoffOD);                    %add min OD value cutoff line
         title('Min OD Data');               %write title
         ylabel('Optical Density');          %write y axis label
         xlabel('Row-Column');               %write x axis label 
      end %function

      function Gh = showGR(obj, curvedata)
      %plot bar graph of GR values for curve data from RowCol table
         set(gcf, 'NextPlot', 'new');        %next plot goes in new figure
         RCtab = curvedata.RowColTable;      %get GR data all rows/cols
         RCtab = sortrows(RCtab,1);          %sort by row values
         initvec = RCtab.GR;                 %extract GR data
         Xlabelvec = obj.getLabelVec([RCtab.ROW, RCtab.COL]); %Xaxis labels 
         Gh = bar(Xlabelvec, initvec);       %bar graph of GR with labels
         title('Growth Rate Data');          %write title
         ylabel('Growth Rate (mu/h)');       %write y axis label
         xlabel('Row-Column');               %write x axis label
      end %function

      function Lh = showLag(obj,curvedata)
      %plot bar graph of lag values from row col table 
         set(gcf, 'NextPlot', 'new');        %next plot goes in new figure
         RCtab = curvedata.RowColTable;      %get lag values from RowCol
         RCtab = sortrows(RCtab,1);          %sort by row vals          
         initvec = RCtab.LAG;                %get lag vector
         Xlabelvec = obj.getLabelVec([RCtab.ROW, RCtab.COL]); %Xaxis labels 
         Lh = bar(Xlabelvec, initvec);            %bar graph of Lag with labels
         title('Lag Time Data');             %write title
         ylabel('Time (h)');                 %write y axis label
         xlabel('Row-Column');               %write x axis lable
      end

      function Mh = showMaxOD(obj, curvedata)
      %plot bar graph of maxOD values as real numbers from row col table
         set(gcf, 'NextPlot', 'new');        %next plot goes in new figure
         RCtab = curvedata.RowColTable;      %get maxOD values from RowCol
         RCtab = sortrows(RCtab,1);          %sort by row vals
         initvec = RCtab.MAX;                %get maxOD vector
         Xlabelvec = obj.getLabelVec([RCtab.ROW, RCtab.COL]); %Xaxis labels
         Mh = bar(Xlabelvec, initvec);       %bar graph of maxOD, labels
         title('Max OD Data');               %write title
         ylabel('Optical Density');          %write y axis label
         xlabel('Row-Column');               %write x axis label
      end %function

      function Th = showPlate(obj)
      %make a tile graph with all 96 growth curves
         close all
         fig = figure(); 
         fig.Position(3:4) =[500, 400];
         t = tiledlayout(8,12,'TileSpacing','none'); %8x12 grid, no spacing   
         xlabel(t,"Time (h)");               %x axis label
         ylabel(t,"Optical Density (600 nm)"); %y axis label
         Ydata = obj.MP_XYY(:,2:end);        %just get Y data values
         ymax_val = max(Ydata,[],"all");     %max Y val for Y axis scale
         for i=2:97                          %for each data col
            nexttile                         %add new tile to graph
            plot(obj.MP_XYY(:,1),obj.MP_XYY(:,i),'.k'); %plot in tile
            ylim([0,ymax_val]);              %set Y limit
            set(gca,'xtick',[],'ytick',[]);  %add tick marks
         end %for
         %assing row and column values
         rowvals = ['A','B','C','D','E','F','G','H'];
         colvals = ["1","2","3","4","5","6","7","8","9","10","11","12"];
         t.Children = flipud(t.Children);    %graph tiles
         for j=1:8                           %label Y axis of tile plot
            rnum = tilenum(t,j,1);           %assign each row value
            ylabel(t.Children(rnum),rowvals(j),'Rotation',0, ...
               'FontWeight','bold');         %set label
         end %for
         for k=1:12                          %X axis values for tile plot
            cnum = tilenum(t,1,k);           %assign column numbers
            title(t.Children(cnum),colvals(k),'FontWeight','bold');
         end %for
         title(t,obj.dataname,'FontSize',12,'FontWeight','bold'); %title
         Th = t;
      end %function

   end %public methods 

   methods (Access = private)
   %-----------------------------------------------------------------------
      function modelXY = getGompertz(~,GrowthRate,Lagtime,MaxOD, ...
            InitVal, TimeVec)
      %From Zwittering et al 1992 paper, Gompertz model from params
      %  Returns: XY matrix (Nx2) with X = time, Y = LnOD for model
         Ndp = length(TimeVec);              %number of timepoints
         modelXY = zeros(Ndp,2);             %results matrix Nx2
         modelXY(:,1) = TimeVec;             %set time values
         start = InitVal;                    %get the initial value
         A = MaxOD + abs(start); %because Aval is negative (Ln) number?
         e = exp(1);                         %assign value of e
         Lambda = Lagtime;                   %assign lag
         mu = GrowthRate;                    %assign growth rate
         modelXY(1,2) = start;               %assign initial Y value
         for i=2:Ndp                         %for each timepoint
            modelXY(i,2) = start + ...       %add predicted val to start pt
               A*exp(-1*exp((mu*e/A)*(Lambda - modelXY(i,1))+1));
         end %for 
      end %function

      function GRdata = getGRdata(obj, rowcolmat, processbyN)
      %Return structure with growth rate params and stats
         LnXYY = obj.getLnData(rowcolmat,obj.blankData.mean); %get LnXYY
         [Ndp, Ncols] = size(LnXYY);         %number of data points, cols
         LnXY = obj.XYY2XY(LnXYY);           %convert XYY to single XY
         Nys = Ncols - 1;                    %number of curves
         istart = 1;                         %initial start index
         iend = processbyN * Nys;            %initial end index
         LR = obj.getLRegdata(LnXY,istart,iend); %initial regression
         startindex = istart;                %save initial start index
         istart = istart + Nys;              %advance start
         iend = iend + Nys;                  %advance end
         while iend <= Ndp                   %check remaining sets  
            LRtemp = obj.getLRegdata(LnXY,istart,iend); %calc GR data
            if LRtemp.Slope > LR.Slope       %check for largest slope
               startindex = istart;          %save new start index
               LR = LRtemp;                  %save regression data
            end
            istart = istart + Nys;           %advance start
            iend = iend + Nys;               %advance end
         end %for
         minOD = obj.getMinOD(LnXY);         %get min LnOD of firstN
         maxOD = obj.getMaxOD(LnXY);         %max LnOD  
         DblTime = 0;                        %preset doubling time
         LagTime = 0;                        %preset lag time
         if LR.Slope > 0                     %check for positive slope
            DblTime = log(2)/LR.Slope;       %calculate doubling time
            if minOD > log(obj.MachineRes)
               LagTime = (minOD - LR.Intercept)/LR.Slope; %calc lag time
            else
               LagTime = 0;      %no lag if min OD <= machine resolution
            end
         end %if 
         %set return strut, get XY model from Gompertz parameters
         modelXY = obj.getGompertz(LR.Slope,LagTime,maxOD,minOD,LnXY(:,1));
         %fill GRdata struct
         GRdata.ObsLnXY = LnXY;              %original XY data Nx2
         GRdata.PredLnXY = modelXY;          %XY data for Gompertz model             
         DataVec = [LR.Slope,DblTime,LagTime,minOD,maxOD]; %GR data vector
         GRdata.ResTable = array2table(DataVec); %make table
         GRdata.ResTable.Properties.VariableNames = {'GR','DBL','LAG',...
            'MIN','MAX'};                    %variable names
         GRdata.RowColMat = rowcolmat;       %row(s) and col(s) used
         GRdata.Npts = Ndp;                  %number of row values
         GRdata.IndexStart = startindex;     %start index for pts used
         GRdata.IndexEnd = startindex + processbyN -1; %end index pts used
         GRdata.Slope = LR.Slope;            %slope of the regression line
         GRdata.YIntercept = LR.Intercept;   %Y intercept regression line
         GRdata.Rsq = LR.Rsquared;           %regression Rsq
         GRdata.Residuals = LR.Residuals;    %residuals from regression
      end %function

      function StatsTable = getStatsTable(obj,restable)
      %calculate stat values for reps from ResTable generated by
      %  getResultsTable
         GRmean = mean(restable.GR);         %growth rate mean
         GRstd = std(restable.GR);           %growth rate std
         GR = [GRmean, GRstd];               %growth rate vector
         DBLmean = mean(restable.DBL);       %doubling time mean
         DBLstd = std(restable.DBL);         %doubling time std
         DBL = [DBLmean, DBLstd];            %doubling time vector
         LAGmean = mean(restable.LAG);       %lag time mean
         LAGstd = std(restable.LAG);         %lag time std
         LAG = [LAGmean, LAGstd];            %lag time vector
         MinODmean = mean(restable.MIN);     %min optical density (ln) mean
         MinODstdev = std(restable.MIN);     %min optical density (ln) std
         MIN = [MinODmean, MinODstdev];      %min optical density (ln) vec
         MaxODmean = mean(restable.MAX);     %max optical density (ln) mean
         MaxODstdev = std(restable.MAX);     %max optical density (ln) std
         MAX = [MaxODmean, MaxODstdev];      %max optical density (ln) vec
         Xlabelvec = obj.getLabelVec([restable.ROW, restable.COL]); %labels 
         labelstr = join(Xlabelvec,',',1);     %make a string from labels
         LabTable = array2table(labelstr,'VariableNames',{'LABEL'}); %name
         stattab = [MIN, GR, DBL, LAG, MAX]; %combined data matrix (2x5) 
         StatsTable = array2table(stattab);  %convert data to table
         StatsTable.Properties.VariableNames = {'MIN','MINstd','GR',...
            'GRstd','DBL','DBLstd','LAG','LAGstd','MAX','MAXstd'}; %names
         StatsTable = [LabTable, StatsTable]; %combine tables (2 x 6 table)
      end %function

      function ResTable = getResultsTable(obj,rowcolmat, processbyN)
      %get a table of results for each row in rowcolmat
         TempTable = table;                  %make a temporary table                
         Nrc = size(rowcolmat,1);            %get number of curves
         for i=1:Nrc                         %for each row in rowcolmat
            temp = obj.getGRdata(rowcolmat(i,:), processbyN); %get results
            TempTable = [TempTable; temp.ResTable]; %#ok<AGROW> 
         end %for
         RCtab = array2table(rowcolmat,"VariableNames",...
            {'ROW','COL'});                  %make rowcol table
         ResTable = [RCtab, TempTable];      %add to data
         ResTable = sortrows(ResTable,'ROW');
      end %function    

      function XYY = getSubXYY(obj,rowcolmat)
      %get sub matrix by row col location in 3 x 12 microplate
         Nrc = size(rowcolmat,1);            %get number of curves
         XYY = zeros(obj.TPs,Nrc+1);         %preallocate return matrix
         XYY(:,1) = obj.MP_XYY(:,1);         %set X (time) values
         for i = 1:Nrc                       %for each row colmod val
            XYtemp = obj.getXY(rowcolmat(i,1),rowcolmat(i,2)); %get curve
            XYY(:,i+1) = XYtemp(:,2);        %assign Y vals
         end %for
      end %function

      function LnXYY = getLnData(obj,rowcolmat,blankval)
      %subtract blankval from all Y vals then convert to natural log
      %NOTE: remove TPs w/ negative or zero Y values prior to Ln transform
      %  This makes sense because all ODs would be positive values
         XYY = obj.getSubXYY(rowcolmat);
         XYY(:,2:end) = XYY(:,2:end) - blankval; %subtrack blank from Ys
         %replace any value less than machine resolution
         LnXYY = max(XYY,obj.MachineRes);    %replace if below resolution
         %remove rows with a negative Y
         % index = all(XYY(:,2:end)>0,2);    %get row indices with pos Y
         % LnXYY = XYY(index,:);             %keep only positive Y vals 
         LnXYY(:,2:end) = log(LnXYY(:,2:end)); %Ln transform of Y vals
      end %function
      
      function LRegdata = getLRegdata(~,LnXY,starttime,endtime)
      %do linear regression for sub set of data (start to end X vals)
      % using Ln OD values
         dataXvals = LnXY(starttime:endtime,1);  %x vals
         dataYvals = LnXY(starttime:endtime,2);  %use Ln vals
         datamat = [dataXvals dataYvals];    %make Nx2 matrix
         LRegdata = LinReg(datamat);         %get regression data
         %The function LinReg returns a data struct with:
         % stats.Time              1 x n array of time values
         % stats.LogCFUpermL       1 x n array of CFU/ml (survivors)
         % stats.PredCFUpermL      1 x n array of predicted CFU/ml values
         % stats.Slope             slope of the regression line
         % stats.StdErrSlope       standard error for slope estimate
         % stats.Intercept         intercept (initial CFU/ml)
         % stats.StdErrIntercept   std. error for intercept estimate
         % stats.Residuals         error values for each data point
         % stats.Rsquared          regression coefficient (R squared value)
      end %function

      function minOD = getMinOD(obj,XY)
      %get only min of firstN data values in matrix 
         if obj.FirstN == 0                      %don't use firstN
            minOD = min(XY(:,2));            %get min of Y vals
         else
            minOD = min(XY(1:obj.FirstN,2));     %use firstN vals only
         end %if
      end %function

      function maxOD = getMaxOD(~,XY)         
      %get max of all OD values
         maxOD = max(XY(:,2));               %assign max OD value           
      end %function

      function XY = getXY(obj,row,col)
      %Get curve from a perticular row and col in the 8 x 12 microplate
         yvals = squeeze(obj.MP_XYN(row,col,1:end)); %OD vals from row col
         xvals = obj.MP_XYY(:,1);            %get X vals
         XY = [xvals yvals];                 %assemble matrix (Nx2)
      end %function

      function XYall = XYY2XY(~,XYY)
         %return Nx2 with one X col and one Y col with all Y reps 
         [Nrows, Ncols] = size(XYY);         %get size of XYY
         NYcols = Ncols - 1;
         XYlen = Nrows*NYcols;               %XYlen = datapoints * Y cols
         resindex = 1;
         XYall = zeros(XYlen,2);             %result matrix with all zeros
         for i=1:NYcols                      %for each Y col
            for j=1:Nrows                    %each row in a Y col
               XYall(resindex,1) = XYY(j,1);
               XYall(resindex,2) = XYY(j,i+1); %get Y data for row,col
               resindex = resindex + 1;      %advance XYall row index
            end %for
         end %for
         XYall = sortrows(XYall,1);
      end %function

      function setXYN(obj, blankval)
         %make XYN from XYY data, subtract blank value from Y data
         colindex = 2;                       %temp col index, first Y col
         for i=1:8                           %each row
            for j=1:12                       %each col 
               %subtract blank value form XYY OD data for each Y col
               % and load XYN matrix
               obj.MP_XYN(i,j,:) = obj.MP_XYY(:,colindex)-blankval; %load Y
               colindex = colindex + 1;      %advance to next col (2-97)
            end %for
         end %for 
      end %function

      function rowcolmat = getRowColMat(~,rowvec, colvec)
      %  convert vector of rows + vector of cols to a row-col matrix (Nx2)
      %  whichever is longer is length of Nx2
      %  if unequal length only the first value of the shorter vec is used
         NRvals = length(rowvec);                  %get length row vector
         NCvals = length(colvec);                  %get length col vector
         Nvals = NRvals * NCvals;                  %number of return rows
         rowcolmat = zeros(Nvals,2);               %return matrix size
         Rindex = 1;                               %temp index for rows
         Cindex = 1;                               %temp index for cols
         for i = 1:Nvals                           %for each entry
            rowcolmat(i,:) = [rowvec(Rindex), colvec(Cindex)]; %make vector
            Rindex = Rindex + 1;                   %advance row index
            if Rindex > NRvals                     %if exceed row max:
               Rindex = 1;                         %reset row index
               Cindex = Cindex + 1;                %advance col index
            end %if
         end %for
      end %function

      function closegraph(~)
      %close all graphs
         close all;                          %close all graphs
      end

      function Xlab = getLabelVec(obj,rowcolvector)
         Xlabels = num2str(rowcolvector);
         Nvals = size(Xlabels,1);
         for i=1:Nvals
            Xlabels(i,1) = obj.setRowLetter(Xlabels(i,1)); 
         end
         Xlab = erase(string(Xlabels), " ");
      end
         
      function Xchar = setRowLetter(~,charval)   
         switch charval
            case '1'
               Xchar = 'A';
            case '2'
               Xchar = 'B';
            case '3'
               Xchar = 'C';
            case '4'
               Xchar = 'D';
            case '5'
               Xchar = 'E';
            case '6'
               Xchar = 'F';
            case '7'
               Xchar = 'G';
            case '8'
               Xchar = 'H';
         end
      end %function

   end %methods -----------------------------------------------

end %class def
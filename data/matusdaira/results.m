% plots conditional means for Matsudaira(2008)
% E[Y|X=x,W=x] using local linear regression
% Y is outcome variable
% W is treatment status {0,1}
% X is forcing variable

%this version: 06/03/2014
%previous versions: 05/21/2014, 04/07/2014


% Y is outcome variable
% W is treatment status {0,1}
% X is the negative of the forcing variable minus cutoff

% VARLIST
% mgrade01 - grade level for math in 01
% ssatyn01 - attended summer school 01
% mdcut01 - math score rel to cutoff in
% rdcut01 - reading score rel to cutoff in
% mret02 =1 if held back (retained) between 01 and 02
% zmscr02 - 2002 math z-score
% zrscr02 - 2002 reading z-score
% mmscr00 - 2000 math score (missing values imputed)
% mrscr00 - 2000 reading score (missing values imputed)


clear all
pathpaper='';


%% file name prefix for saving output
filenamepre='results_matsu';



%%
load grade5
n=length(mdcut01);


%enter each outcome variable you wanna try as different column of YY
%Y=rdge1-rdge0;
%YY=[zrscr02,zmscr02];
YY=zmscr02;
Yname={'math'};

%enter each forcing variable you wanna try as different column of WW
XX=min([rdcut01 mdcut01],[],2);
Xname={'min(reading,math)'};

%to fit the standard RD framework, eligibility for treatment occurs
%when X>=0, and not X<=0; this will matter for RDD3
XX=-XX;
W=ssatyn01;



%figure counter
fig=1;

%number of grid points for plotting
np=[100 100]; %LHS and RHS of cutoff
      

%table that keeps jumpt estimates, s.e. and t-stats
%for conditional mean, conditional mean given treatment, and probability
%jump
%for all cobinations of outcome and forcing variables
table=zeros((size(XX,2)*3),size(YY,2),4);


% %optimal bandwidth (Imbens and Kalyanaraman)
% [est,se,h_all]=rd_optbandwidth(Y,X,W,zeros(n,1),0,0,0);

for y=1:size(YY,2);
    for x=1:size(XX,2);
        Y=YY(:,y);
        X=XX(:,x);


        %% ************************************************************************
        %ploting the conditional mean E[Y|X]
        %**************************************************************************




        %[late,se,h,bias]=rd_optbandwidth(Y,X,W,zeros(n,1),0,0,0);
        [jump,se,~,h,extras]=RDD(Y,X,(X>=0),0);

%         gridL=linspace(min(X),0,np(1));
%         gridR=linspace(0,max(X),np(2));
        
        gridL=linspace(-70,0,np(1));
        gridR=linspace(0,70,np(2));

        %removes those points from the grid for which 
        %there won't be enough obs to run the LLR with bandwidth h
        nL=length(X(X<0));
        nR=length(X(X>0));
        gridL=gridL(sum(abs((repmat(X(X<0),1,np(1))-repmat(gridL,nL,1))/h)<1)>=10);
        gridR=gridR(sum(abs((repmat(X(X>0),1,np(2))-repmat(gridR,nR,1))/h)<1)>=10);


        figure(fig);
        fig=fig+1;
        plot(gridL, LLR(gridL',h,Y(X<0),X(X<0)),'blue',gridR,LLR(gridR',h,Y(X>0),X(X>0)),'red');
        title(['E(Y|X), h=',num2str(h)])
        
        ylabel(['outcome var: ',Yname{y}]);
        %xlabel(['forcing var: ',Xname{x}, '; jump:', num2str(jump), '; t-stat: ',num2str(jump/se),'; p-value: ',num2str(normcdf(-abs(jump/se))*2)]);
        xlabel(['forcing var: ',Xname{x}, ';   jump:', num2str(jump), ';   S.E.: ',num2str(se)]);
        
        %cropping margins of the figure
        tightInset = get(gca, 'TightInset'); %TightInset defines the tightest bounding box that encloses the axes and its labels and title.
        position(1) = tightInset(1);
        position(2) = tightInset(2);
        position(3) = 1 - tightInset(1) - tightInset(3);
        position(4) = 1 - tightInset(2) - tightInset(4);
        set(gca, 'Position', position);
        
        %positioning the figure for printing
        set(gcf,'Units','inches');
        position=get(gcf,'Position');
        set(gcf, 'PaperPosition', [0 0 position(3) position(4)]); %This places plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [position(3) position(4)]); %Set the paper to have width 5 and height 5.
        
        %saves in eps and pdf formats
        saveas(gcf, [pathpaper,filenamepre,'_EY_X_case',num2str(y),'_',num2str(x)],'eps');
        saveas(gcf, [pathpaper,filenamepre,'_EY_X_case',num2str(y),'_',num2str(x)],'pdf');







        %% ************************************************************************
        %ploting the conditional mean E[Y|X,W=w] for w=0,1
        %**************************************************************************
        estw=zeros(2,1); sew=estw; hw=estw; biasw=estw;

        for w=0:1


            %[estw(w+1),sew(w+1),hw(w+1),biasw(w+1)]=rd_optbandwidth2(Y(W==w),X(W==w),(X(W==w)<=0),zeros(sum(W==w),1),0,0,0);
            %[estw(w+1),sew(w+1),biasw(w+1),hw(w+1)]=RDD2(Y,X,W,0,w);
            [estw(w+1),sew(w+1),biasw(w+1)]=RDD2(Y,X,W,0,w,h);


            %number of grid points for plotting
%             gridL=linspace(min(X(W==w)),0,np(1));
%             gridR=linspace(0,max(X(W==w)),np(2));

            gridL=linspace(-70,0,np(1));
            gridR=linspace(0,70,np(2));

            %removes those points from the grid for which 
            %there won't be enough obs to run the LLR with bandwidth h
            nL=length(X((X<0)&(W==w)));
            nR=length(X((X>0)&(W==w)));
            %gridL=gridL(sum(abs((repmat(X((X<0)&(W==w)),1,np(1))-repmat(gridL,nL,1))/hw(w+1))<1)>=10);
            %gridR=gridR(sum(abs((repmat(X((X>0)&(W==w)),1,np(2))-repmat(gridR,nR,1))/hw(w+1))<1)>=10);
            gridL=gridL(sum(abs((repmat(X((X<0)&(W==w)),1,np(1))-repmat(gridL,nL,1))/h)<1)>=10);
            gridR=gridR(sum(abs((repmat(X((X>0)&(W==w)),1,np(2))-repmat(gridR,nR,1))/h)<1)>=10);

            figure(fig);
            fig=fig+1;
            %plot(gridL, LLR(gridL',hw(w+1),Y((X<0)&(W==w)),X((X<0)&(W==w))),'blue',gridR,LLR(gridR',hw(w+1),Y(X>0&W==w),X(X>0&W==w)),'red');
            plot(gridL, LLR(gridL',h,Y((X<0)&(W==w)),X((X<0)&(W==w))),'blue',gridR,LLR(gridR',h,Y(X>0&W==w),X(X>0&W==w)),'red');
            %title(['E(Y|W=',num2str(w),',X), h=',num2str(hw(w+1))])
            title(['E(Y|W=',num2str(w),',X), h=',num2str(h)])
            ylabel(['outcome var: ',Yname{y}]);
            %xlabel(['jump: ',num2str(estw(w+1)),'; t-stat: ',num2str((estw(w+1)-biasw(w+1))/sew(w+1)),'; p-value: ',num2str(normcdf(-abs((estw(w+1)-biasw(w+1))/sew(w+1)))*2)]);
            xlabel(['forcing var: ',Xname{x}, ';   jump: ',num2str(estw(w+1)),';   S.E.: ',num2str(sew(w+1))]);
            %print(figh,'-depsc',['EY_X_W',num2str(w),'_case',num2str(y),'_',num2str(x),'.eps']);
            
            %cropping margins of the figure
            tightInset = get(gca, 'TightInset'); %TightInset defines the tightest bounding box that encloses the axes and its labels and title.
            position(1) = tightInset(1);
            position(2) = tightInset(2);
            position(3) = 1 - tightInset(1) - tightInset(3);
            position(4) = 1 - tightInset(2) - tightInset(4);
            set(gca, 'Position', position);

            %positioning the figure for printing
            set(gcf,'Units','inches');
            position=get(gcf,'Position');
            set(gcf, 'PaperPosition', [0 0 position(3) position(4)]); %This places plot at left hand corner with width 5 and height 5.
            set(gcf, 'PaperSize', [position(3) position(4)]); %Set the paper to have width 5 and height 5.

            %saves in eps and pdf formats
            saveas(gcf, [pathpaper,filenamepre,'_EY_X_W',num2str(w),'_case',num2str(y),'_',num2str(x)],'eps');            
            saveas(gcf, [pathpaper,filenamepre,'_EY_X_W',num2str(w),'_case',num2str(y),'_',num2str(x)],'pdf');            

        end

        %% ************************************************************************
        %ploting the conditional probability prob(W=1|X)   
        %**************************************************************************


        %[estp,sep,biasp,hp]=RDD(W,X,(X>=0),0);
        [estp,sep,biasp]=RDD(W,X,(X>=0),0,h);




        %gridL=linspace(min(X),0,np(1));
        %gridR=linspace(0,max(X),np(2));
        gridL=linspace(-70,0,np(1));
        gridR=linspace(0,70,np(2));

        %removes those points from the grid for which 
        %there won't be enough obs to run the LLR with bandwidth h
        nL=length(X(X<0));
        nR=length(X(X>0));
        %gridL=gridL(sum(abs((repmat(X(X<0),1,np(1))-repmat(gridL,nL,1))/hp)<1)>=10);
        %gridR=gridR(sum(abs((repmat(X(X>0),1,np(2))-repmat(gridR,nR,1))/hp)<1)>=10);
        gridL=gridL(sum(abs((repmat(X(X<0),1,np(1))-repmat(gridL,nL,1))/h)<1)>=10);
        gridR=gridR(sum(abs((repmat(X(X>0),1,np(2))-repmat(gridR,nR,1))/h)<1)>=10);

        figure(fig);
        fig=fig+1;
        %plot(gridL, LLR(gridL',hp,W(X<0),X(X<0)),'blue',gridR,LLR(gridR',hp,W(X>0),X(X>0)),'red');
        plot(gridL, LLR(gridL',h,W(X<0),X(X<0)),'blue',gridR,LLR(gridR',h,W(X>0),X(X>0)),'red');
        %title(['P(W=1|X), h=',num2str(hp)]);
        title(['P(W=1|X), h=',num2str(h)]);
        %xlabel(['jump: ',num2str(estp),'; t-stat: ',num2str((estp-biasp)/sep),'; p-value: ',num2str(normcdf(-abs((estp-biasp)/sep))*2)]);
        %xlabel(['forcing var: ',Xname{x}, '; jump: ',num2str(estp),'; t-stat: ',num2str(estp/sep),'; p-value: ',num2str(normcdf(-abs(estp/sep))*2)]);
        xlabel(['forcing var: ',Xname{x}, ';   jump: ',num2str(estp),';   S.E.: ',num2str(sep)]);
        ylim([-.0001 1.0001]);
        %print(figh,'-depsc',['EW_X_case',num2str(y),'_',num2str(x),'.eps']);

        %cropping margins of the figure
        tightInset = get(gca, 'TightInset'); %TightInset defines the tightest bounding box that encloses the axes and its labels and title.
        position(1) = tightInset(1);
        position(2) = tightInset(2);
        position(3) = 1 - tightInset(1) - tightInset(3);
        position(4) = 1 - tightInset(2) - tightInset(4);
        set(gca, 'Position', position);

        %positioning the figure for printing
        set(gcf,'Units','inches');
        position=get(gcf,'Position');
        set(gcf, 'PaperPosition', [0 0 position(3) position(4)]); %This places plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [position(3) position(4)]); %Set the paper to have width 5 and height 5.

        %saves in eps and pdf formats
        saveas(gcf, [pathpaper,filenamepre,'_EW_X_case',num2str(y),'_',num2str(x)],'eps');            
        saveas(gcf, [pathpaper,filenamepre,'_EW_X_case',num2str(y),'_',num2str(x)],'pdf');  


        % %% ************************************************************************
        % %ploting the density f(X)
        % %**************************************************************************
        % 
        % %kernel function
        % ker = @(x) (abs(x)<=1).*(1-abs(x));
        % 
        % oh=2; %factor that multiplies Silverman's bandwidth
        % 
        % 
        % %bandwidth - Silverman's rule for Edge Kernel and Normal X  
        % hs=oh*2.57*std(X)*(n^(-1/5));
        % 
        % gridL=linspace(min(X),max(X(X<0)),np(1));
        % gridR=linspace(min(X(X>0)),max(X),np(2));
        % 
        % 
        % 
        % %removes those points from the grid for which 
        % %there are not enough obs around 
        % nL=length(X(X<0));
        % nR=length(X(X>0));
        % gridL=gridL(sum(abs((repmat(X(X<0),1,np(1))-repmat(gridL,nL,1))/hs)<1)>=10);
        % gridR=gridR(sum(abs((repmat(X(X>0),1,np(2))-repmat(gridR,nR,1))/hs)<1)>=10);
        % 
        % np2=[length(gridL),length(gridR)];
        % 
        % figh=figure(fig);
        % fig=fig+1;
        % %you need to multiply the density estimator by 2 because you're using
        % %X>0 or X<0 obs only
        % plot(gridL, (2/(n*hs))*sum(ker( (repmat(X(X<0),1,np2(1))-repmat(gridL,nL,1))  /hs  )),'blue',...
        %      gridR, (2/(n*hs))*sum(ker( (repmat(X(X>0),1,np2(2))-repmat(gridR,nR,1))  /hs  )),'red');
        % title(['f(X), h=',num2str(hs)])
        % print(figh,'-depsc','f_X.eps');
        %         
        %  
        % 
        % 
        % 
        % 
        % 
        % %% ************************************************************************
        % %ploting the conditional densities f(X|W)
        % %**************************************************************************
        % 
        % %bandwidths
        % hsw=zeros(2,1);
        % 
        % for w=0:1 
        % 
        %     
        %     nw=length(X(W==w));
        % 
        %     %bandwidth - Silverman's rule for Edge Kernel and Normal X  
        %     
        % 
        %     hsw(w+1)=oh*2.57*std(X(W==w))*(nw^(-1/5));
        %     gridL=linspace(min(X(W==w)),max(X((W==w)&(X<0))),np(1));
        %     gridR=linspace(min(X((X>0)&(W==w))),max(X(W==w)),np(2));
        % 
        %     
        %     
        %     %removes those points from the grid for which 
        %     %there are not enough obs around 
        %     nL=length(X((X<0)&(W==w)));
        %     nR=length(X((X>0)&(W==w)));
        %     gridL=gridL(sum(abs((repmat(X((X<0)&(W==w)),1,np(1))-repmat(gridL,nL,1))/hsw(w+1))<1)>=10);
        %     gridR=gridR(sum(abs((repmat(X((X>0)&(W==w)),1,np(2))-repmat(gridR,nR,1))/hsw(w+1))<1)>=10);
        % 
        %     np2=[length(gridL),length(gridR)];
        %     
        %     figh=figure(fig);
        %     fig=fig+1;
        %     %you need to multiply the density estimator by 2 because you're using
        %     %X>0 or X<0 obs only
        %     plot(gridL, (2/(nw*hsw(w+1)))*sum(ker(((repmat(X((W==w)&(X<0)),1,np2(1)))-repmat(gridL,nL,1))/hsw(w+1))),'blue',...
        %          gridR, (2/(nw*hsw(w+1)))*sum(ker(((repmat(X((W==w)&(X>0)),1,np2(2)))-repmat(gridR,nR,1))/hsw(w+1))),'red');
        %     title(['f(X|W=',num2str(w),'), h=',num2str(hsw(w+1))])
        %     print(figh,'-depsc',['f_X_W',num2str(w),'.eps']);
        %         
        %  
        % end
        % 

        %% ************************************************************************
        %table of results: jump in conditional means E[Y|X] and E[W|X] 
        %**************************************************************************
        fprintf('-----------------------------------------------------------\n');
        fprintf(['Outcome var: ',Yname{y},'; forcing var: ',Xname{x},'\n']);
        fprintf('Bandwidth: %-9.4g\n',h); 


        fprintf('-----------------------------------------------------------\n');

        fprintf('E[W|X=0+]-E[W|X=0-]: %-9.4g \n',estp);
        fprintf('s.e. : %-9.4g \n',sep);
        fprintf('t-stat : %-9.4g \n',estp/sep);

        fprintf('-----------------------------------------------------------\n');


        fprintf('E[Y|X=0+]-E[Y|X=0-]: %-9.4g \n',jump);
        fprintf('s.e. : %-9.4g \n',se);
        fprintf('t-stat : %-9.4g \n',jump/se);






        table([1;2;3]+(x-1)*3,y,1)=[jump;se;jump/se];
        table([1;2;3]+(x-1)*3,y,4)=[estp;sep;estp/sep];

        %% ************************************************************************
        %table of results: jump in conditional means E[Y|X,W=w] and E[W|X,W=w] 
        %**************************************************************************

        fprintf('-----------------------------------------------------------\n');

        for w=0:1

            fprintf(['E[Y|X=0+,W=',num2str(w),']-E[Y|X=0-,W=',num2str(w), ']: %-9.4g \n'],estw(w+1));
            fprintf('s.e. : %-9.4g \n',sew(w+1));
            fprintf('t-stat : %-9.4g \n',estw(w+1)/sew(w+1));

            fprintf('-----------------------------------------------------------\n');
            table([1;2;3]+(x-1)*3,y,2+w)=[estw(w+1);sew(w+1);estw(w+1)/sew(w+1)];

        end




        %% ************************************************************************
        %table of results: conditional means by compliance group
        %**************************************************************************
        [ eta, cov, bias] = RDD3 (Y,X,W,0,h);
        %eta - (4x1) estimated
        %       E[Y(0)|W=c,complier];
        %       E[Y(0)|W=c,nevertaker];
        %       E[Y(1)|W=c,complier];
        %       E[Y(1)|W=c,alwaystaker];
        names={'E[Y(0)|X=0,complier]','E[Y(0)|X=0,nevertaker]','E[Y(1)|X=0,complier]','E[Y(1)|X=0,alwaystaker]'};

        
        fprintf([names{3},' - ', names{1},': %-9.4g \n'],eta(3)-eta(1));
        fprintf('s.e. : %-9.4g \n',sqrt(cov(3,3)+cov(1,1)-2*cov(1,3)));
        fprintf('t-stat : %-9.4g \n',(eta(3)-eta(1))/sqrt(cov(3,3)+cov(1,1)-2*cov(1,3)));

        fprintf('-----------------------------------------------------------\n');
        fprintf([names{4},' - ', names{3},': %-9.4g \n'],eta(4)-eta(3));
        fprintf('s.e. : %-9.4g \n',sqrt(cov(4,4)+cov(3,3)-2*cov(3,4)));
        fprintf('t-stat : %-9.4g \n',(eta(4)-eta(3))/sqrt(cov(4,4)+cov(3,3)-2*cov(3,4)));

        fprintf('-----------------------------------------------------------\n');        
        
        fprintf([names{2},' - ', names{1},': %-9.4g \n'],eta(2)-eta(1));
        fprintf('s.e. : %-9.4g \n',sqrt(cov(2,2)+cov(1,1)-2*cov(1,2)));
        fprintf('t-stat : %-9.4g \n',(eta(2)-eta(1))/sqrt(cov(2,2)+cov(1,1)-2*cov(1,2)));

        fprintf('-----------------------------------------------------------\n');
        
        

        for c=1:2
            fprintf([names{c},': %-9.4g \n'],eta(c));
            fprintf('s.e. : %-9.4g \n',sqrt(cov(c,c)));
            fprintf('t-stat : %-9.4g \n',eta(c)/sqrt(cov(c,c)));

            fprintf('-----------------------------------------------------------\n');

        end



        for c=3:4
            fprintf([names{c},': %-9.4g \n'],eta(c));
            fprintf('s.e. : %-9.4g \n',sqrt(cov(c,c)));
            fprintf('t-stat : %-9.4g \n',eta(c)/sqrt(cov(c,c)));

            fprintf('-----------------------------------------------------------\n');

        end



    end
end



%**************************************************************************
%Export results into Latex files 
%**************************************************************************


%**************************************************************************
%- filenamepre_table_case.tex -
%one table for each case of forcing and outcome
%variable combination; this table includes includes:  
%jump in conditional probability of treatment status,
%jump in conditional mean of Y, 
%jump in conditional mean of Y given treatment status;
%- in each one, we report the jump and se for all combinations
%of forcing and outcome variables;
%**************************************************************************


    


fid = fopen([pathpaper,filenamepre,'_table_case',num2str(y),'_',num2str(x),'.tex'], 'w');




    
fprintf(fid,'\\begin{table}[H]\n');
fprintf(fid,'\\begin{center}\n');


fprintf(fid,'\\caption{Matsudaira Data - Optimal Bandwidth $h_{opt}=%-9.4g$}\n',h);
fprintf(fid,'\\label{tab:matsu_table}\n');


fprintf(fid, '\\begin{tabular}{c|ccc}\n');
fprintf(fid, 'Outcome & Estimand & Estimate & S.E. \\\\ \n');
fprintf(fid, '\\hline  \n');
fprintf(fid, 'Summer School & $\\mathbb{E}[W_i(1)-W_i(0)|X_i=0]$ & %-9.4g & (%-9.4g) \\\\ \n',[estp;sep]);
fprintf(fid, 'Math Score & $\\mathbb{E}[Y_i(1)-Y_i(0)|X_i=0]$ & %-9.4g & (%-9.4g) \\\\ \n',[jump;se]);
fprintf(fid, 'Math Score & $\\mathbb{E}[Y_i(1)-Y_i(0)|X_i=0,G_i=c]$ & %-9.4g & (%-9.4g) \\\\ \n',[eta(3)-eta(1);sqrt(cov(3,3)+cov(1,1)-2*cov(1,3))]);
fprintf(fid, 'Math Score & $\\mathbb{E}[Y_i(1)|X_i=0,G_i=a] - \\mathbb{E}[Y_i(1)|X_i=0,G_i=c]$ & %-9.4g & (%-9.4g) \\\\ \n',[eta(4)-eta(3);sqrt(cov(4,4)+cov(3,3)-2*cov(3,4))]);
fprintf(fid, 'Math Score & $\\mathbb{E}[Y_i(0)|X_i=0,G_i=n] - \\mathbb{E}[Y_i(0)|X_i=0,G_i=c]$ & %-9.4g & (%-9.4g) \\\\ \n',[eta(2)-eta(1);sqrt(cov(2,2)+cov(1,1)-2*cov(1,2))]);

fprintf(fid, '\\end{tabular}\n');
fprintf(fid,'\\end{center}\n');

fprintf(fid,'\\end{table}\n\n\n');
fclose(fid);




%**************************************************************************
%- filenamepre_tables.tex -
%one table for each: jump in conditional mean of Y, 
%conditional mean of Y given treatment status, and conditional
%probability;
%- in each one, we report the jump and se for all combinations
%of forcing and outcome variables in XX and YY;
%**************************************************************************



fid = fopen([pathpaper,filenamepre,'_tables.tex'], 'w');

%tables 1 to 4


for tab=1:4
    %eval(['table',num2str(tab),'=table(:,:,',num2str(tab),');']);
    
    fprintf(fid,'\\begin{table}[H]\n');
    fprintf(fid,'\\begin{center}\n');

    if tab==1
        fprintf(fid,'\\caption{Jump in Conditional Mean}\n');
    elseif (tab==2)||(tab==3)
        fprintf(fid,['\\caption{Jump in Conditional Mean Given Treatment=',num2str(tab-2),'}\n']);


    elseif tab==4
        fprintf(fid,'\\caption{Jump in Conditional Probability}\n');
        
    end
    fprintf(fid,'\\label{table%d}\n',tab);

    fprintf(fid, '\\begin{tabular}{cc|');
    for y=1:size(YY,2)
        fprintf(fid,'c');
    end
    fprintf(fid,'}\n');
    fprintf(fid,' & &\\multicolumn{%d}{c}{Outcome Variable}\\\\\n',size(YY,2));
    fprintf(fid,'Forcing var. & ');
    for y=1:size(YY,2)
        fprintf(fid,[' &',Yname{y}],y);
    end
    fprintf(fid,'  \\\\\n');
    fprintf(fid,'  \\hline \n');

    for x=1:size(XX,2);
        fprintf(fid,' & jump ');
        for y=1:size(YY,2);
            fprintf(fid,' & %-9.4g ',table(1+(x-1)*3,y,tab));
        end
        fprintf(fid,'  \\\\\n');
        fprintf(fid,[Xname{x},' & se ']);
        for y=1:size(YY,2);
            fprintf(fid,' & %-9.4g ',table(2+(x-1)*3,y,tab));
        end
        fprintf(fid,'  \\\\\n');        

        fprintf(fid,' & t-stat ');
        for y=1:size(YY,2);
            fprintf(fid,' & %-9.4g ',table(3+(x-1)*3,y,tab));
        end
        fprintf(fid,'  \\\\\n');
        fprintf(fid,'  \\hline \n');

    end
    

    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid,'\\end{center}\n');

    fprintf(fid,'\\end{table}\n\n\n');

end


fclose(fid);

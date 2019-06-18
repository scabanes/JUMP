%Programma per creare la griglia polare sia in px che in cm %%% Input: un fes%% Output: - [ Grid_Xp, Grid_Yp]  punti della griglia polare in px con              %% Nraggi e Ncerchi definiti dall'utente centrata in X0=0,Y0=0               %% cio� l'origine della figura (estremo sinistro basso), mi              %% serve per l'interpolazione dei fes   %%      - [Grid_Xp_cm, Grid_Yp_cm] punti della griglia polare in cm   %%        centrata in r=0, theta=0     %%      - r_max_cm, raggio massimo in cm    %%      - vettore del raggio da 0 a r_max_cm con passo dr fid2=fopen([roots,Name,'/fes2_.FES'])%[a,data]=readFesInfo('fes2_esp1.FES') LEGGE HEADER    data1=fread(fid2, 1, 'int');    disp(sprintf('headerLenght %f',data1))    %%pause        data2=fread(fid2, 1, 'int');    disp(sprintf('maxFeatures %f',data2))        data3=fread(fid2, 1, 'ushort');    disp(sprintf('ImgRows %f',data3))        data4=fread(fid2, 1, 'ushort');    disp(sprintf('ImgCols %f',data4))        data5=fread(fid2, 1, 'ushort');    disp(sprintf('RoiX0 %f',data5))        data6=fread(fid2, 1, 'ushort');    disp(sprintf('RoiY0 %f',data6))          data7=fread(fid2, 1, 'ushort');    disp(sprintf('RoiWidth %f',data7))        data8=fread(fid2, 1, 'ushort');    disp(sprintf('RoiHeight %f',data8))          data9=fread(fid2, 1, 'ushort');        disp(sprintf('StartFr %f',data9))      %pause    data10=fread(fid2, 1, 'ushort');    disp(sprintf('Step %f',data10))     %pause        data11=fread(fid2, 1, 'int');    disp(sprintf('FrameNumber %f',data11))    %pause    fseek(fid2,-data1,'cof');    fseek(fid2,data1,'bof');X0 = data7/2;Y0 = data8/2;% Ncerchi_in = input('Number of circles (small=40 o big=60): ','s');% Ncerchi = sscanf(Ncerchi_in,'%f');Ncerchi=Nri;% Nraggi_in = input('Number of rays (es. 120): ','s');% Nraggi = sscanf(Nraggi_in,'%f');Nraggi=Nti;% cmpx_in = input(' Inserire il cmpx dell esp : ','s');% cmpx = sscanf(cmpx_in,'%f');cmpx=0.066; MaxRadius = X0; %%� in pixel           Rpoints = 0:MaxRadius/Ncerchi:MaxRadius;    Rpoints = Rpoints(1:end-1);    Tpoints = 0:2*pi/Nraggi:2*pi;    Tpoints = Tpoints(1:end-1);        [GridR,GridT]=meshgrid(Rpoints,Tpoints);        %% punti della griglia polare in px con centro in r=0, theta=0    Grid_Xp = GridR.*cos(GridT);    Grid_Yp = GridR.*sin(GridT);        %%% punti della griglia polare in px con centro in X0=0,Y0=0, estremo    %%% sinistro basso    %%% mi serve per l'interpolazione dei fes        Grid_Xp = Grid_Xp + X0+data5;    Grid_Yp = Grid_Yp + Y0+data6;            %% trasformazione dei punti della griglia in cm con origine della    %% griglia in r=0, theta=0% alcuni esp in vasca piccola%cmpx =0.045;  Grid_Xp_cm = (Grid_Xp-X0-data5)*cmpx;Grid_Yp_cm = (Grid_Yp-Y0-data6)*cmpx;    r_max_cm = max(Grid_Xp_cm(:));r = Rpoints*cmpx;dr = (MaxRadius/Nri)*cmpx;dtheta = (2.*pi)/Nti;theta = [0:dtheta:2.*pi-dtheta];% Beta-effectf0=2.*Omega;s = Omega.^2/(2.*g);Hc = H0-(s*(Lx.*Lx + Ly.*Ly)/12);H=Hc+(r.^2).*(s);
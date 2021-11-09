%% Interpolation and resampling of raw kinematic data - The output of this set of calculations is saved in Trajectories.mat
% 'Raw_Trajectories'  is a matlab struct containing all raw kinematic data
% for a particular subject.'
% 'Corrected onset' contains the manually-selected timing of movement onset for the respective trial used to crop data during post-processing.
% 'Corrected offset'contains the manually-selected timing of movement offset for the respective trial used to crop data during post-processing.
% 'TNumber'  is an integer referencing the index of the trial to be processed.

% Function description
% Performs processing steps on kinematic data for a single trial as described in Furmanek et al.
% Authors: Mariusz P. Furmanek, Madhur Mangalam, Kyle Lockwood, Mathew Yarossi, and Gene Tunik
% Contact:  m.furmanek@northeastern.edu
% Date:     May 2021
% Version:  1

% Extract trial data
Input_Data = Raw_Trajectories{TNumber};

% Interpolation parameters
fs = 100;   % 100 Hz sampling frequency

% Filtering parameters
fc = 6;         % Cutoff frequency
ord = 4;        % Order
type = 'low';   % Type

Start = find(Input_Data{1}==Corrected_Onset(TNumber)) - 37; % Subtract 37 to include the trajectory from 500 ms before movement onset 
End = find(Input_Data{1}==Corrected_Offset(TNumber));
Original_Time = Input_Data{1}(Start:End)-Input_Data{1}(Start);

for Chan = 2:length(Input_Data) % Iterate over all channels in raw trial data
    
    Data_to_Interp = Input_Data{Chan}(Start:End)*100; % Convert m to cm
    
    % Interpolation: 75 Hz raw data is upsampled to 100 Hz
    Interp_Time = [1/fs:1/fs:Original_Time(end)/1000]';
    Interp_Data = interp1(Original_Time/1000,Data_to_Interp,Interp_Time);
    
    % Filtering: apply Butterworth filter to raw data
    [b,a] = butter(ord,fc/(fs/2),type);
    Output_Data(:,Chan) = filtfilt(b,a,Interp_Data);
    
end

Output_Data(:,1) = [-490:10:Interp_Time(end - 50)*1000];

Resampled = Output_Data;

%% Kinematic profiles - The output of this set of calculations is saved in Profiles.mat

% Aperture, Aperture velocity and Aperture acceleration

Thumb_X = Output_Data(50:end,5); Thumb_Y = Output_Data(50:end,6);
Finger_X = Output_Data(50:end,8); Finger_Y = Output_Data(50:end,9);
Aper_Tmp = sqrt((Finger_X - Thumb_X).^2 + (Finger_Y - Thumb_Y).^2);
Aper = Aper_Tmp(2:end); % So that aperture and aperture velocity have the same number of samples

AperVelo_Tmp = diff(Aper_Tmp)*100;
fs = 100; [b,a] = butter(4,6/(fs/2)); % Apply 4th-order 6 Hz low-pass filter
AperVelo = filtfilt(b,a,AperVelo_Tmp);

AperAcc_Tmp = diff([0;AperVelo])*100;
fs = 100; [b,a] = butter(4,6/(fs/2)); % Apply 4th-order 6 Hz low-pass filter
AperAcc = filtfilt(b,a,AperAcc_Tmp);

% Tranport distance, Transport velocity and Transport acceleration

Wrist_X = Output_Data(50:end,2); Wrist_Y = Output_Data(50:end,3);
Wrist_X0 = Output_Data(50,2); Wrist_Y0 = Output_Data(50,3);
TransDis_Tmp = sqrt((Wrist_X - Wrist_X0).^2 + (Wrist_Y - Wrist_Y0).^2);
TransDis = TransDis_Tmp(2:end);

TransVelo_Tmp = diff(TransDis_Tmp)*100;
fs = 100; [b,a] = butter(4,6/(fs/2)); % Apply 4th-order 6 Hz low-pass filter
TransVelo = filtfilt(b,a,TransVelo_Tmp);

TransAcc_Tmp = diff([0;TransVelo])*100;
fs = 100; [b,a] = butter(4,6/(fs/2)); % Apply 4th-order 6 Hz low-pass filter
TransAcc = filtfilt(b,a,TransAcc_Tmp);

%% Kinematic features - The output of this set of calculations is saved in Features.mat

% Movement time, MT [ms]
MT = sum(~isnan(TransDis),1)*10; K = MT;

% Peak transport velocity, Peak_TV [cm/s] and
[Peak_TV T_Peak_TV_Tmp] = max(TransVelo);

% Time to peak transport velocity, T_Peak_TV [ms]
T_Peak_TV = T_Peak_TV_Tmp*10;

% Peak transport acceleration, Peak_TA [cm/s2]
[Peak_TA T_Peak_TA_Tmp] = max(TransAcc);
Peak_TA =  Peak_TA;

% Time to peak transport velocity, T_Peak_TV [ms]
T_Peak_TA = T_Peak_TA_Tmp*10;

% Peak transport decelersation, Peak_TD [cm/s2]
[Peak_TD T_Peak_TD_Tmp] = min(TransAcc);

% Time to peak transport deceleration, T_Peak_TD [cm/s2]
T_Peak_TD = T_Peak_TD_Tmp*10;

% Peak aperture, Peak_A [cm]
[Peak_A OT_Tmp] = max(Aper);

% Peak aperture velocity, Peak_AV [cm/s]
[Peak_AV T_Peak_AV_Tmp] = max(AperVelo);

% Time to peak aperture velocity, T_Peak_AV [ms]
T_Peak_AV = T_Peak_AV_Tmp*10;

% Peak aperture acceleration, Peak_AA [cm/s2]
[Peak_AD T_Peak_AD_Tmp] = min(AperAcc);

% Time to peak aperture acceleration, T_Peak_AA [ms]
T_Peak_AD = T_Peak_AD_Tmp*10;

% Opening time, OT [ms]
OT = OT_Tmp*10;

% Closing time, CT [ms]
CT = (sum(~isnan(TransDis),1) - OT_Tmp)*10;

% Opening distance, OD [cm]
Ind = OT_Tmp + [0:20:20*(size(TransDis,2) - 1)];
OD = TransDis(Ind);

% Closure distance, CD [cm]
Ind_2 = sum(~isnan(TransDis),1) + [0:200:200*(size(TransDis,2) - 1)];
CD = TransDis(Ind_2) - TransDis(Ind);

% Transport velocity at closure onset, TV_CO [cm/s]
TV_CO = TransVelo(Ind);

% Transport acceleration at closure onset, TA_CO [cm/s2]
TA_CO = TransAcc(Ind);
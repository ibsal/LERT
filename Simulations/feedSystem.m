function out = feedSystem(in)
%   in.




end

%% useful subfunctions
ipaKinV(Temperature)
NitrousVaporPressure(Temperature)

%% Models Needed:
n2oKinV(Temperature) % For calculation of RE
pidPDrop(in.ox.pid, Mdot?) % This will need to be coupled with the engine model, initial velocity would also be weird as an input?




Feed system states: 

Nitrous temperature
Nitrous pressure 
Nitrous vapor pressure 
IPA temperature
IPA pressure 
MdotNitrous
MdotIPA

% Fluid line object structure. Starts from tank section, ends at inlet to
% engine 
in.ox.pid = [["LINE", length, roughness, diameter], ["COMPONENT", diameter, CD]];
in.fuel.pid = [["LINE", length, roughness, diameter], ["COMPONENT", diameter, CD]];
in.vessel.oxTemp
in.vessel.ipaTemp
in.geom.oxHeight
in.geom.fuelHeight
in.geom.id


out.geom.oxHeight
out.geom.fuelHeight

% Fill process
% Install fixed volume of IPA at stp
IPAVolume = ID^2 * 0.25 * pi * IPAHeight;
% Flip so IPA is on top
% Begin fill of saturated N2O 
N2OVolume = ID^2 * 0.25 * pi * (TotalHeight-IPAHeight);
N2OPressure
% Fill terminates once remaining volume is availabel 
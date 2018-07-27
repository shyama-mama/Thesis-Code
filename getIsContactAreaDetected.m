%% Helper function
% Author: Heba Khamis
% Date: 06/10/2015
%%

function isContactAreaDetected = getIsContactAreaDetected(contactAreaFrame)

temp = find(contactAreaFrame==1,1,'first');
isContactAreaDetected = ~isempty(temp);

end
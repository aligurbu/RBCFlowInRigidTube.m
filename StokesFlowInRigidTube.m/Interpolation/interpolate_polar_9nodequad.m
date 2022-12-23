function Ntilde = interpolate_polar_9nodequad(rho, thet, s01, s02)
%%
% Computes (N(s01+rho*cos(theta),s02+rho*sin(theta)) - N(s01,s02))/rho
% Compute interpolation functions in polar coordinates  
% for 9-node quadrilateral element
%%
ct = cos(thet);
st = sin(thet);

rst = rho*st;
rct = rho*ct;
rctst = rho*ct*st;

S01_p_s01_st = s01*(s01+1)*st;
S01_m_s01_st = s01*(s01-1)*st;
S02_p_s02_ct = s02*(s02+1)*ct;
S02_m_s02_ct = s02*(s02-1)*ct;

S01_st = (1-s01^2)*st;
S02_ct = (1-s02^2)*ct;

ss01 = 2*s01;
ss02 = 2*s02;

s01_p = ss01+1;
s01_m = ss01-1;
s02_p = ss02+1;
s02_m = ss02-1;

s01_m_rct = s01_m + rct;
s02_m_rst = s02_m + rst;

rctst_s01_m_rct = rctst *(s01_m_rct);
rctst_s01_p_p_rct = rctst *(s01_p + rct);
rctst_ss01_p_rct = rctst *(ss01 + rct);

Ntilde = [( S01_m_s01_st*(s02_m_rst) ...
          + S02_m_s02_ct*(s01_m_rct) ...
          + rctst_s01_m_rct*(s02_m_rst))/4;  
          ( S01_p_s01_st*(s02_m_rst) ...
          + S02_m_s02_ct*(s01_p + rct) ...
          + rctst_s01_p_p_rct*(s02_m_rst))/4;  
          ( S01_p_s01_st*(s02_p + rst) ...
          + S02_p_s02_ct*(s01_p + rct) ...
          + rctst_s01_p_p_rct*(s02_p + rst))/4;
          ( S01_m_s01_st*(s02_p + rst) ...
          + S02_p_s02_ct*(s01_m_rct) ...
          + rctst_s01_m_rct*(s02_p + rst))/4;
          ( S01_st*(s02_m_rst) ...
          - S02_m_s02_ct*(ss01 + rct) ...
          - rctst_ss01_p_rct*(s02_m_rst))/2;
          (-S01_p_s01_st*(ss02 + rst) ...
          + S02_ct*(s01_p + rct) ...
          - rctst_s01_p_p_rct*(ss02 + rst))/2;
          ( S01_st*(s02_p + rst) ...
          - S02_p_s02_ct*(ss01 + rct) ...
          - rctst_ss01_p_rct*(s02_p + rst))/2;
          (-S01_m_s01_st*(ss02 + rst) ...
          + S02_ct*(s01_m_rct) ...
          - rctst_s01_m_rct*(ss02 + rst))/2;
          - S01_st*(ss02+rst) ...
          - S02_ct*(ss01+rct) ...
          + rctst_ss01_p_rct*(ss02+rst)];

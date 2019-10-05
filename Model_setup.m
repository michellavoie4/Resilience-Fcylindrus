%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curation of the F cylindrus metabolic network
% based on the model of P. tricornutum of Levering et al. (2016) and Broddrick
% et al. (2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;  %  Close previous program windows.
clear all;  %  Clear previous program variables.
clc;  %  Clear the screen of previous text.
format short; % Easier to see numbers

% time the code
tic;

% Initiate Cobra Toolbox if needed
% initCobraToolbox()

changeCobraSolver('gurobi', 'all')
% Assert that gurobi is installed and working properly
if ans == 0
error('Error loading gurobi')
end

% Do not read The Initial Network in iLB1025.xml as it introduces error in the
% cobra model. e.g., No subsystems are present
% model = readCbModel('iLB1025.xml')

% You can also save and load the model in a .mat file 
% 1) save the model in a .mat file
% save Fragila.mat model
% 2) load the .mat file
% load('Fragila.mat', 'model')

% Either, load directly the model of F. cylindrus as a .mat file from
% levering et al (2016)
load('iLB1025bon.mat')

% Since there are errors when loading the model with readCbModel()
% We must run this function. Otherwise the function getModelSubSystems does
% not work.
pti = convertOldStyleModel(pti)
model = pti
model= generateRules(model)

model.id = 'Adapted from iLB1025'
model.description = 'Metabolic network reconstruction of Fragilariopsis cylindrus'

% Remove the field 'metChEBIID', which does not match the required
% properties of the function verifyModel()
model = rmfield(model, 'metChEBIID')

% this work, the model can be verified
verifyModel(model, 'simpleCheck', true)
% writeCbModel(pti, 'fileName', 'modelptitest.xml')

% Curation of 'model'
mets_length = length(model.mets)
rxns_length = length(model.rxns)

% Updating the grRules field with the genes of F. cylindrus
% model.grRules = {}
% Fixing the bugs related to reactions deletion : When deleting reactions,
% grRules becomes 'double' rather than 'character'
model.grRules()
GPR_length = length(model.rxns) % gr.Rules must be of the same length than model.rxns
cellarray = cell(GPR_length,1)
% convert a cell array (double) to character
cellArrayChar = cellfun(@num2str, cellarray, 'UniformOutput', false)
% change model1.grRules to character
model.grRules = cellArrayChar 

iscellstr(model.grRules) % True if it is a cell array of character array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding all gene numbers or protein IDs of Fragilariopsis cylindrus CCMP1102, which is
% the same strain than CCMP3323.

model.grRules{4} = '(207264)'  % Do not add {}, otherwise it is no longer a cell arrays of character array, but it becomes a cell array with one cell array in each cell
model.grRules{310} = '(277092)'
%model.grRules{1768} = No genes, partial protein
model.grRules{1} = '(242337 and 210953 and 208175 and 268240 and 267742)'
model.grRules{381} = '(185720 and 180463)'
model.grRules{2} = '(219585)'
model.grRules{3} = '(219585)'
model.grRules{5} = '(163305) or (143564) or (153439)'
model.grRules{6} = '(259666)'
model.grRules{311} = '(207456) or (260715)'
model.grRules{7} = '(158325)' % Manually annotated
model.grRules{274} = '(277211)' % Manually annotated
model.grRules{8} = '(229452) or (208444)' 
model.grRules{9} = '(239808)' % Manually annotated
% model.grRules{1910} = No genes targeted to the chloroplast
model.grRules{10} = '(226043)'
model.grRules{11} = '(170472) or (209612)'
model.grRules{12} = '(238714)' % Manually annotated
model.grRules{14} = '(193396)'
model.grRules{13} = '(291486)'
model.grRules{15} = '(206834) or (208444)'
model.grRules{16} = '(194532)'
% model.grRules{17} = No genes targeted to mitochondria
model.grRules{315} = '(183228) or (188279)'
model.grRules{18} = '(179924)'
model.grRules{1883} = '(223701)'
model.grRules{19} = '(209415)'
model.grRules{20} = '(275301) or (245629)'  % Manually annotated
model.grRules{21} = '(291486)'
model.grRules{22} = '(291486)'
model.grRules{23} = '(277202)'  % Manually annotated
model.grRules{24} = '(274454) or (235074)'  % Manually annotated
model.grRules{25} = '(291634)'
model.grRules{26} = '(178758)'
model.grRules{27} = '(206834) or (208444)'
model.grRules{28} = '(275554)'
model.grRules{29} = '(187273) or (268822)'
model.grRules{30} = '(187750)'
model.grRules{31} = '(226762)'
model.grRules{32} = '(212647)'
model.grRules([33,34,35],1) = {'(291486)'}
model.grRules{36} = '(224559)'
model.grRules{37} = '(207939)'
model.grRules{475} = '(180647)'
model.grRules{1145} = '(276098)'
model.grRules{38} = '(269372)'  % Manually annotated
model.grRules{316} = '(274452)' 
model.grRules{386} = '(275634)'  % Manually annotated 
model.grRules{39} = '(268008 and 260490 and 235520)'
model.grRules{40} = '(186177)'
model.grRules{41} = '(263231)'
model.grRules{42} = '(183259)'  % Manually annotated
model.grRules{405} = '(186955)'  % Manually annotated
model.grRules([43,44,45,46],1) = {'(291486)'}
model.grRules{47} = '(201716)'
model.grRules{48} = '(239548) or (252263)'
model.grRules{314} = '(291505) or (277738)'
model.grRules{49} = '(241035) or (229638)'  % Manually annotated
model.grRules{53} = '(259279)'
model.grRules{50} = '(274784)' % Manually annotated
model.grRules{55} = '(242838) or (224993) or (260072)'  % Manually annotated, I added one more mitochondrial PK, 3 instead of 2
model.grRules{389} = '(228971)'  % Manually annotated
model.grRules{51} = '(206063) or (274784) or (225883)'  % Manually annotated
model.grRules{52} = '(206063) or (274784) or (225883)'  % Manually annotated
model.grRules{54} = '(225939) or (247496)'  % Manually annotated
model.grRules{56} = '(247304)'
% model.grRules{57} = No chloroplastic analogous genes
model.grRules{58} = '(268008 and 260490 and 235520)'
model.grRules{59} = '(206063) or (274784) or (225883)'  % Manually annotated
model.grRules{60} = '(157897)'
model.grRules{61} = '(269865) or (157897)' 
model.grRules{62} = '(211817)'
% model.grRules{63} = No significant similarities.
model.grRules{64} = '(206834) or (208444)' 
model.grRules{65} = '(186500)'
model.grRules{66} = '(271375)'  % Manually annotated
model.grRules{73} = '(263399)'
model.grRules{74} = '(187383)'  % Manually annotated
model.grRules{197} = '(263399)'
model.grRules{198} = '(263399)'
model.grRules{199} = '(180836)'
model.grRules{354} = '(180836)'
model.grRules{67} = '(275493)'
model.grRules{68} = '(213551)'
model.grRules{69} = '(185594) or (188245)'
model.grRules{70} = '(275358)'
model.grRules{71} = '(275358)'
model.grRules{72} = '(262390)'
% model.grRules{75} = 'No mitochondrial genes, as rxn 17
model.grRules{76} = '(277191) or (268209)'   % Manually annotated
model.grRules{77} = '(246065)'
model.grRules{423} = '(236986)'
% model.grRules{78} = No chloroplastic genes
model.grRules{79} = '(268008 and 260490 and 235520)'
model.grRules{80} = '(232330)'
model.grRules{1568} = '(180557)'
model.grRules{1571} = '(247580)'
model.grRules{81} = '(157897)'
model.grRules{82} = '(183955)'
model.grRules{83} = '(183955)'
model.grRules{84} = '(209911) or (177704)'   % Manually annotated
model.grRules{85} = '(192807)'
model.grRules{365} = '(166658)'
model.grRules{86} = '(170146)'
model.grRules{87} = '(269867)'
model.grRules{88} = '(270428)'   % Manually annotated
model.grRules{90} = '(209715)'   % Manually annotated
model.grRules{89} = '(209933)'
model.grRules{320} = '(209502) or (208673)'   % Manually annotated
model.grRules{387} = '(210020)'   % Manually annotated
model.grRules{91} = '(261136)'
model.grRules{92} = '(184892)'
model.grRules{93} = '(207660)'   % Manually annotated
model.grRules{94} = '(166620)'
% model.grRules{95} = no genes, partial
model.grRules{96} = '(180330)'
model.grRules{97} = '(180330)'
model.grRules{98} = '(259901)'
% model.grRules{99} = no genes, partial
model.grRules{100} = '(205852)'     % Manually annotated
model.grRules([101:104],1) = {'(157897)'} 
% model.grRules{105} = No chloroplastic genes
% model.grRules{106} = No homologous genes
model.grRules{107} = '(192565)'
model.grRules{108} = '(185965)'
model.grRules{109} = '(259321)'
model.grRules{110} = '(267731)'
model.grRules{111} = '(213551)'
model.grRules{112} = '(275493)'
model.grRules{113} = '(271472 and 173934 and 227880 and 196661)'
model.grRules{114} = '(268607)'
model.grRules{115} = '(271881)'  % Manually annotated
model.grRules{116} = '(226968)'
model.grRules{117} = '(202604)'
model.grRules{118} = '(277320)'  % Manually annotated
% model.grRules{119} = No genes, partial sequence
model.grRules{120} = '(242870)'
model.grRules{121} = '(225004)'
model.grRules{122} = '(225004)'
model.grRules([123,124], 1) = {'(211512)'}
model.grRules([125,155:157,294:297,309],1) = {'(173527)'}
model.grRules([590,594,598],1) = {'(186602)'}
model.grRules([1320,1324,1328,1332,1336,1340,1344],1) = {'(232542)'}
model.grRules{126} = '(268473)'
model.grRules{127} = '(268500) or (209393)'   % Manually annotated
model.grRules{128} = '(259321)'  % Manually annotated
model.grRules{129} = '(267731)'
model.grRules{130} = '(174198) or (158118)'   % Manually annotated
model.grRules{131} = '(269285) or (206976)'
model.grRules{236} = '(269285) or (206976)'
model.grRules{338} = '(191311) or (182808)'
model.grRules{132} = '(269005)'
model.grRules([133,134], 1) = {'(287526)'}
% model.grRules{135,136} = No chloroplastic genes
model.grRules([137,138], 1) = {'(263072)'}   % Manually annotated
% model.grRules{139} = No homologous genes
model.grRules([140,217,471:474,491:492,579:582,1346,1993,1994],1) = {'(206834) or (208444)'}
model.grRules([992:997],1) = {'(215629)'}
model.grRules{141} = '(259810)'  % Manually annotated
model.grRules{1315} = '(259810)'
model.grRules{142} = '(276895)'  % Manually annotated
model.grRules{143} = '(206583)' 
model.grRules{144} = '(267781)'  % Manually annotated
model.grRules{145} = '(160102)'
% model.grRules{146} = No genes, partial
model.grRules{147} = '(213551)'
model.grRules{148} = '(275493)'
model.grRules{149} = '(277095)'  % Manually annotated.
model.grRules{150} = '(167330)'
model.grRules{325} = '(207931)'
model.grRules{151} = '(149227)'
model.grRules{152} = '(172194)'
model.grRules{153} = '(225230) or (172086)'
model.grRules{154} = '(276978) or (269015)'
model.grRules{158} = '(273154) or (254523) or (273154)'
model.grRules{329} = '(189638)'
model.grRules{159} = '(259678)'
model.grRules{160} = '(259678)'
model.grRules{161} = '(191403)'
model.grRules{162} = '(267911)'
model.grRules{163} = '(167448)'   % Manually annotated
model.grRules{164} = '(268162)'   % Manually annotated
model.grRules{165} = '(186354 and 208648)' 
model.grRules{166} = '(186354 and 208648)'
model.grRules{174} = '(186354 and 208648)'
model.grRules{167} = '(275954)'
model.grRules{168} = '(291511)'
model.grRules{169} = '(197039) or (197110)'
model.grRules{170} = '(166620)'   % Manually annotated
model.grRules{171} = '(166620)'   % Manually annotated
% model.grRules{172} = No analogous genes
model.grRules{173} = '(202663)'   % Manually annotated
model.grRules([175,176:182],1) = {'(274965) or (260660)'}
model.grRules([587,591,595],1) = {'(273692)'}
model.grRules([1317, 1321, 1325, 1329, 1333, 1337, 1341],1) = {'(228745)'}
model.grRules{189} = '(273718)'
model.grRules{424} = '(211032)'
model.grRules{183} = '(228674)'
model.grRules{184} = '(207678)'
% model.grRules{185} = No chloroplastic genes
% model.grRules{186, 187} = No significant similarities
model.grRules{188} = '(226908)'
% model.grRules{190} = No cytosolic genes
% model.grRules{191} = No genes, partial
model.grRules{192} = '(223582)'
model.grRules{193} = '(223582)'
model.grRules{194} = '(211031)'
model.grRules{195} = '(269031)'
model.grRules{196} = '(228745) or (260660)'
model.grRules{200} = '(207847)'
model.grRules{201} = '(190753)'
model.grRules{202} = '(179359)'
model.grRules{203} = '(186585)'
model.grRules{204} = '(209964)'
model.grRules{205} = '(163305) or (143564) or (153439)'
model.grRules{206} = '(207280)'
model.grRules{207} = '(187013)'
model.grRules{208} = '(238951)'  % Manually annotated
model.grRules{209} = '(179077)'
model.grRules{210} = '(209097)'
model.grRules{211} = '(275387)'
model.grRules{212} = '(246826)'  % Manually annotated
model.grRules{213} = '(187807)'  % Manually annotated
model.grRules{214} = '(255725)'
model.grRules{215} = '(212988)'
model.grRules{216} = '(196819)'
model.grRules{218} = '(261055 and 170289 and 228557)'
model.grRules{219} = '(268444)'  % Manually annotated
% model.grRules{220} = No significat homologous genes
model.grRules{221} = '(189108)'
% model.grRules{222} = No significant homolog
model.grRules{223} = '(267640)'  % Manually annotated
model.grRules{224} = '(140996)'
model.grRules{225} = '(274501)'
model.grRules{226} = '(264173)'  % Manually annotated
model.grRules{227} = '(205101)'
model.grRules{228} = '(218256)'
model.grRules{229} = '(267397) or (240805) or (267395)'
model.grRules{230} = '(261082)'  % Manually annotated
model.grRules{231} = '(209097)'
model.grRules{232} = '(195266)'
model.grRules{233} = '(185542)'
model.grRules{234} = '(261984)'
model.grRules{235} = '(170693)'
model.grRules{237} = '(195411)' 
model.grRules{238} = '(180836)'
model.grRules{239} = '(232580)'  % Manually annotated
model.grRules{240} = '(226008)'  %% ?? boarderline chloroplast
model.grRules{241} = '(170693)'
% model.grRules{242} = No homologous
model.grRules{243} = '(270634)'
model.grRules{244} = '(182520) or (228581)'  % Manually annotated
% model.grRules{1961} = No homologous
model.grRules{245} = '(226164)'  % Manually annotated
model.grRules{246} = '(238951)'
model.grRules{247} = '(260963)'
model.grRules{248} = '(291551)'  % Manually annotated
model.grRules{249} = '(268822)' 
model.grRules{250} = '(224738) or (182516)' 
model.grRules{251} = '(267684)'
model.grRules{252} = '(196266)'
model.grRules{253} = '(186177)'
model.grRules{254} = '(263231)'
model.grRules{255} = '(178356)'
model.grRules{256} = '(208511)'
model.grRules{257} = '(185193) or (206934)'
model.grRules{258} = '(208713) or (274773)'
% model.grRules{1708} = No chloroplastic genes, CYSTL_h
model.grRules{259} = '(267136)'  % Boarderline mitochondrial
model.grRules{287} = '(206884)'
model.grRules{260} = '(224334)'
% model.grRules{1721} = No chloroplastic genes, HISTL_h
model.grRules{261} = '(167622) or (276205)'
model.grRules{262} = '(212866)'
model.grRules{1714} = '(187452)'
% model.grRules{263} = No mitochondrial genes, LYSTL_m
model.grRules{358} = '(168292)'
model.grRules{264} = '(205100)'
% model.grRules{265} = No mitochondrial genes, PHETL_m
model.grRules{1723} = '(197215)'
model.grRules{1724} = '(192847)'
model.grRules{266} = '(275420)'
model.grRules{267} = '(231803)'
model.grRules{268} = '(271163)'
model.grRules{269} = '(268525)'  % Manually annotated, TRPTL_c
% model.grRules{1730} = No chloroplastic genes, TRPTL_h
model.grRules{270} = '(258963)'
% model.grRules{1727} = No mitochondrial genes, TYRTL_m
model.grRules{271} = '(225934)'
% model.grRules{272} = No chloroplastic genes, VALTL_h
model.grRules{273} = '(167852)'
model.grRules{1877} = '(167852)'
model.grRules{274} = '(268562) or (208582) or (188055)' % 268562 is manually annotated
model.grRules{391} = '(209949)'
model.grRules{276} = '(209986)'
model.grRules{277} = '(223701)'
model.grRules{279} = '(173652) or (211467)'
model.grRules{280} = '(211179) or (239120)'
model.grRules{388} = '(238967)'
model.grRules{281} = '(195174)'
model.grRules{282} = '(195174)'
model.grRules{1637} = '(186879) or (205346)'
model.grRules{283} = '(195174)'
model.grRules{284} = '(195174)'
model.grRules{285} = '(267260)'  % Manually annotated
model.grRules{286} = '(250461)'
model.grRules{288} = '(188848 and 278067 and 273772 and 239694 and 259027 and 211896 and 234556 and 191854 and 145259 and 171138 and 170047 and 268581 and 183388 and 211175)' % Complex 1 of mitochondrial electron transport. In reality, this has more than 40 subunits!
model.grRules([289, 291, 467, 468, 2151], 1) = {'(195641)'}
model.grRules{1345} = '(215839)'
model.grRules{290} = '(271375 and 224733 and 275698 and 276117)'  % Three out of 4 are manually annotated
model.grRules{292} = '(194752)'
model.grRules{293} = '(230231)'
model.grRules([298,299:305,308,593,597,1319,1323,1327,1331,1335,1339,1343], 1) = {'(232938) or (273974)'}
model.grRules{589} = '(240917)'
model.grRules{306} = '(259810)'
model.grRules{1316} = '(259810)'
model.grRules{307} = '(274583)'
model.grRules{470} = '(145710)'
% model.grRules{1314} = No gene
model.grRules{312} = '(261588) or (249103)'  % Manually annotated, CAT_1, CAT_2
model.grRules{317} = '(274749)'
model.grRules{318} = '(207790)' 
model.grRules{313} = '(170890)'
model.grRules{319} = '(208253)'
model.grRules{321} = '(196723)'
model.grRules{322} = '(196723)'
model.grRules{1770} = '(249281)'
% model.grRules{1771} = No genes, partial
% model.grRules{323} = No genes, partial
model.grRules{324} = '(206899)'  % Manually annotated DXPS_h
model.grRules{326} = '(185189)'  % Manually annotated GGPS_c
model.grRules{327} = '(268690)'  % Manually annotated UDPSQS_h
% model.grRules{2096} = No genes in P. tricornutum model, so no homologous
model.grRules{328} = '(263086)' 
model.grRules{330} = '(206407)' 
model.grRules{331} = '(206407)' 
model.grRules{332} = '(274061)'
model.grRules{333} = '(274061)'
model.grRules{1756} = '(274061)'
model.grRules{334} = '(179475)'
model.grRules{335} = '(185143)'
model.grRules{336} = '(277974)'
% model.grRules{337} = No mitochondrial genes
% model.grRules([339:347],1) = No chloroplastic genes
model.grRules([588,592,596],1) = {'(208887)'}
model.grRules([1318,1322,1326,1330,1334,1338,1342],1) = {'(197230)'}
model.grRules{348} = '(182884)'
model.grRules{349} = '(196763)'
model.grRules{350} = '(250137) or (249291)'  % Manually annotated, HOXG_c
model.grRules{351} = '(223502)'  % Manually annotated, CHLPAS_h
model.grRules{352} = '(169676)'
model.grRules{372} = '(169676)'
% model.grRules{353} = No genes in P. tricornutum, so no homologous
model.grRules{355} = '(188279)'
model.grRules{356} = '(185157)'
% model.grRules{1703} = No chloroplastic genes found for SERTL_h
model.grRules{357} = '(229124)'
model.grRules{1696} = '(184752)'
model.grRules{359} = '(185131)'
model.grRules{360} = '(234949)'
% model.grRules{361} = No similarities.
% model.grRules{362} = 'No cytosolic enzymes found
model.grRules{364} = '(275453)'
% model.grRules{1593} = No chloroplastic genes found
model.grRules{366} = '(268965)'
% model.grRules{438} = No chloroplastic genes found
model.grRules{367} = '(223273)' % Manually annotated, IMPHL_1_h
model.grRules{368} = '(223273)' % Manually annotated, IMPHL_2_h
model.grRules{369} = '(262171)'
model.grRules{370} = '(259614)' % Manually annotated, HMBDPO_h
model.grRules{371} = '(277748)'
model.grRules{373} = '(273775)'
model.grRules{374} = '(273775)'
%%%%% in 378, 379 and 380 : Check proteins 500.....
model.grRules{378} = '(172395 and 167955 and 206783 and 235080 and 268077 and 212328 and 206118 and 223681)' % One protein manually annotated
model.grRules{380} = '(267711 and 165544 and 193422 and 200350 and 169285 and 177731)' % All manually annotated
model.grRules{379} = '(268792)'
model.grRules{382} = '(243753)' % manually annotated, GLCt_e
model.grRules{383} = '(183952)'
model.grRules{384} = '(216851)'
model.grRules{385} = '(263182) or (273844)'
% model.grRules{390} = No genes, PPTPt_x. This reaction will be potentially deleted
% since it lead to a deadend
model.grRules{392} = '(178753)'
model.grRules{393} = '(209949)'
model.grRules{394} = '(268562)'
% model.grRules{395} = No genes, partial
model.grRules{396} = '(259433)'
% model.grRules{397, 398, 401, 402} = No genes in P. tricornutum
model.grRules{399} = '(178586)' % Manually annotated, EDD_m
model.grRules{400} = '(267632)' % Manually annotated, EDA_m
% Check rxn 403: the numbers 500...
model.grRules([403, 1452], 1) = {'(184116 and (147482 or 196615 or 195439)'}
model.grRules{404} = '(232978)' % Manually annotated, ACK_c
model.grRules{406} = '(168262)' % Manually annotated, IDH_m
model.grRules{407} = '(177813)'
model.grRules{408} = '(182884)'
model.grRules{409} = '(192073)'
model.grRules{410} = '(209960)'
model.grRules{1559} = '(209458)'
% model.grRules{411, 418, 477, 686, 688, 951, 952, 2077, 2113} = No genes in P. tricornutum
model.grRules{412} = '(208511)'
model.grRules{413} = '(268766)'
model.grRules{414} = '(224163)'
model.grRules{415} = '(233094)'
model.grRules{2137} = '(269031)'
model.grRules{416} = '(244936)'
model.grRules{417} = '(235647)' % Manually annotated, UREASE_m
% model.grRules{419} = No cytosolic genes
model.grRules{420} = '(189498) or (225672)'
% model.grRules{421} = No mitochondrial genes
model.grRules{422} = '(195514)' 
model.grRules{425} = '(247064)'
model.grRules{426} = '(273694)'
% model.grRules{427} = No chloroplastic genes
model.grRules{428} = '(207067)'
% model.grRules{429} = No genes, partial
% model.grRules{431} = No chloroplastic genes
model.grRules{432} = '(227748)'
model.grRules{433} = '(259082)'
model.grRules{434} = '(262715)'
model.grRules{435} = '(226930)'
model.grRules{436} = '(168185)'
model.grRules{437} = '(192916)'
model.grRules{438} = '(168185)'
model.grRules{439} = '(234764)'
model.grRules{440} = '(187917)'
model.grRules{441} = '(226194)'
model.grRules{442} = '(177779)'
model.grRules{443} = '(177779)'
model.grRules{444} = '(223475)'
model.grRules{445} = '(178139)'
model.grRules{446} = '(223475)'
% model.grRules{447} = No chloroplastic genes found
% model.grRules{1755} = No mitochondrial genes found
% model.grRules{448} = No chloroplastic genes found
model.grRules{449} = '(259813)'
model.grRules{1453} = '(229511)'
model.grRules{1454} = '(275677)'
model.grRules{450} = '(259813)'
model.grRules{1455} = '(229511)'
model.grRules{1456} = '(275677)'
model.grRules{278} = '(259813)'
model.grRules{1457} = '(229511)'
model.grRules{1458} = '(275677)'
model.grRules{451} = '(180245)' % Manually annotated, DHSR_h
% model.grRules{452} = No genes, partial
model.grRules{453} = '(238714) or (245452)' % ADDITION: 245452 is Manually annotated, PPHTA_h 
model.grRules{454} = '(191710)'  % Manually annotated, ARODH_h
model.grRules{455} = '(188526)'  % Manually annotated, PHEH_c
model.grRules{456} = '(212618)'
model.grRules{457} = '(225343)'
model.grRules{458} = '(180417)'
model.grRules{459} = '(179168)'
model.grRules{460} = '(178087)'
model.grRules{461} = '(179413)'
% model.grRules{1670} = No chloroplastic genes found
model.grRules{462} = '(143661)'
model.grRules{463} = '(238411)'
model.grRules{464} = '(233147)'
% model.grRules{469} = No genes, partial
model.grRules{1915} = '(233147) or (274583)'
model.grRules{465} = '(169535)'
model.grRules{466} = '(169535)'
model.grRules{476} = '(180508)'
model.grRules{478} = '(277748)'
model.grRules{479} = '(277748)'
model.grRules{480} = '(277748)'
model.grRules{481} = '(190693)'
model.grRules{482} = '(223160)'
model.grRules{483} = '(264100)'
model.grRules{484} = '(264100)'
model.grRules{485} = '(209086)'  % Manually annotated, G1PUT_c
model.grRules{486} = '(237242)'
model.grRules{487} = '(213392)'  % Manually annotated, UGP_PGM_c
model.grRules{488} = '(222153) or (141920)' 
model.grRules{489} = '(271986)' 
model.grRules{490} = '(270903)'
model.grRules([493:496, 2148], 1) = {'(180963)'}
% model.grRules{915,916,917,1347,2000} = No genes, partial
model.grRules([497,498:524,2097,2149], 1) = {'(205304)'}
% 918 to 924, 2081, 2088 = No genes, partial
model.grRules{1348} = '(189436)'
model.grRules([525:548,641:649,2098], 1) = {'(256838)'}
model.grRules([925,926,927,928,929,930,2082,2089,2094], 1) = {'(208572)'}
model.grRules([549:560,2150], 1) = {'(180315)'}
model.grRules([561,562], 1) = {'(277802)'}
model.grRules([563,564,567,568], 1) = {'(228533)'}
model.grRules([607,608,609,610], 1) = {'(268672) or (226788)'}
% 618, 619: No genes, partial
% 565, 566, 569, 570, 605, 606, 611, 612 : No genes, partial
% 572 to 578 , 1842 : No genes in P. tricornutum
model.grRules([986:991], 1) = {'(250000)'}  % Manually annotated
model.grRules([583:586,613,614,696:713,903:913], 1) = {'(182179) or (247727) or (205304)'} 
model.grRules{599} = '(177742)' 
model.grRules{600} = '(177742)'
model.grRules{615} = '(182629) or (277401)'
model.grRules([601:604,616,617,714:734,2102], 1) = {'(208289)'} 
model.grRules([620:640,694], 1) = {'(247182)'} 
model.grRules{650} = '(229525)'
model.grRules([946,947], 1) = {'(264457)'} 
model.grRules{1349} = '(229525)'
model.grRules([651:683], 1) = {'(178092)'} 
model.grRules{684} = '(180315)'
model.grRules{685} = '(212958)'
model.grRules{687} = '(180508)'
model.grRules{689} = '(226382)'
model.grRules{690} = '(210853)'
model.grRules{691} = '(182698)'
model.grRules{2153} = '(194512)'
model.grRules{2154} = '(203452) or (209729)'
model.grRules{2152} = '(167010) or (178564)'
model.grRules{692} = '(248393)'
% model.grRules{693} = No similarities
model.grRules([736:739], 1) = {'(195685)'} 
% model.grRules{740, 741} = No genes, partial
model.grRules([742:770,2106,2107], 1) = {'(276226) or (181385) or (182857)'} 
model.grRules([975:980], 1) = {'(276226)'} 
model.grRules([771:835,2108:2112], 1) = {'(189188)'} 
model.grRules([836:864,865:897,898:902,1976:1983,2114:2116], 1) = {'(248141) or (154996) or (182179) or (240289) or (246374) or (205304)'} 
model.grRules([903:913], 1) = {'(229252) or (247727)'} 
model.grRules{931} = '(275222)' % Manually annotated, GALE_h 
model.grRules{985} = '(170931)' % Manually annotated,  GALE_c
model.grRules([932:934, 2090], 1) = {'(275442)'} 
model.grRules([935, 936, 2085, 2092], 1) = {'(178606)'} 
model.grRules([937, 938, 2086], 1) = {'(225029)'} 
model.grRules{939} = '(188442)'
model.grRules([964:973, 2084, 2087, 2091, 2093], 1) = {'(185101)'} 
model.grRules{2083} = '(275442)'
model.grRules([940,941], 1) = {'(194825)'} 
model.grRules([942:945,2095], 1) = {'(207999)'} 
model.grRules{2099} = '(262455)' 
model.grRules{2100} = '(178606)' 
model.grRules{948} = '(244906)' 
% model.grRules{949, 950} = No genes, partial
model.grRules{1350} = '(185515)' 
model.grRules{953} = '(177929)'
model.grRules([954:959], 1) = {'(138778)'} 
model.grRules([960:963], 1) = {'(232054)'} 
% model.grRules{983} = No genes, partial
model.grRules{984} = '(209655)'  % Manually annotated, G3MGGH_c
model.grRules([998,1002,1006,1010,1014,1018,1022,1026,1030,1038,1044,1052,1055,1063,1069,1077,1083], 1) = {'(210789)'} % Manually annotated, ACOAOX240
model.grRules([999,1003,1007,1011,1015,1019,1023,1027,1031,1035,1041,1045,1049,1056,1060,1066,1070,1074,1080,1086], 1) = {'(207194)'} % Manually annotated, ECH_HADH 
model.grRules([1147,1151,1155,1159,1163,1167,1171,1175,1179,1184,1188,1192,1197,1201,1205,1209,1216,1223,1227,1231,1235,1241,1245,1249,1253,1258,1262,1266,1270,1277,1284,1290,1294,1298,1304,1308], 1) = {'(183437) or (202663)'}  %  Both manually annotated
model.grRules([1000, 1004, 1008, 1012, 1016, 1020, 1024, 1028, 1032, 1036, 1042, 1046, 1050, 1057, 1061, 1067, 1071, 1075, 1081, 1087], 1) = {'(207194)'} % Manually annotated, ECH_HADH
model.grRules([1148, 1152, 1156,  1160, 1164, 1168, 1172, 1176, 1180, 1185, 1189, 1193, 1198, 1202, 1206, 1210, 1217, 1224, 1228, 1232, 1236, 1242, 1246, 1250, 1254, 1259, 1263, 1267, 1271, 1278, 1285, 1291, 1295, 1299, 1305, 1309], 1) = {'(183437) or (270026)'} % Both manually annotated
model.grRules([1001, 1005, 1009, 1013, 1017, 1021, 1025, 1029, 1033, 1037, 1043, 1047, 1051, 1058, 1062, 1068, 1072, 1076, 1082, 1088], 1) = {'(207546)'} 
model.grRules([1149, 1153, 1157, 1161, 1165, 1169, 1173, 1181, 1186, 1190, 1194, 1199, 1203, 1207, 1211, 1218, 1225, 1229, 1233, 1237, 1243, 1247, 1251, 1255, 1260, 1264, 1268, 1272, 1279, 1286, 1292, 1296, 1300, 1306, 1310], 1) = {'(207731) or (268162)'}  % 268162 is manually annotated
model.grRules([1034, 1040, 1048, 1054, 1059, 1065, 1073, 1079, 1085], 1) = {'(207194)'} % Manually annotated, ECOA
model.grRules([1182, 1195, 1208, 1214, 1221, 1234, 1240, 1256, 1269, 1275, 1282, 1289, 1297, 1303], 1) = {'(202663)'}  % Manually annotated, Enoyl-Coa hydratase
model.grRules([1039, 1053, 1064, 1078, 1084], 1) = {'(212175)'} 
model.grRules([1213, 1220, 1239, 1274, 1281, 1288, 1302], 1) = {'(175014)'} 
% 1089 : Spontaneous, No genes
model.grRules{1090} = '(179105)'  
model.grRules{1897} = '(156300)'
model.grRules{1091} = '(179105)'
model.grRules{1092} = '(179105)'
model.grRules([1106,1107:1116,1995,1996], 1) = {'(156300)'} 
% model.grRules{1093, 1103, 1104, 1105, 1117 to 1130, 1997, 1998, 2012, 2013} = No genes in P. tricornutum.  
model.grRules{1095} = '(262191)'  % Manually annotated, TMLOX_m
model.grRules{1096} = '(143661)' 
model.grRules{1097} = '(247568)' 
model.grRules{1098} = '(269865)'
model.grRules{1099} = '(179475)'
model.grRules{1100} = '(236986)'
model.grRules{1101} = '(208511)'
model.grRules{1102} = '(234666)' % Manually annotated, GBBOX_m
model.grRules([1131, 1132:1144, 1999,  2000], 1) = {'(179105)'} 
model.grRules([1146, 1150, 1154, 1158, 1162, 1166, 1170, 1178, 1183, 1187, 1191, 1196, 1200, 1204, 1212, 1215, 1219, 1222, 1226, 1230, 1238, 1244, 1248, 1252, 1257, 1261, 1265, 1273, 1276, 1280, 1283, 1287, 1293, 1301, 1307], 1) = {'(212170) or (209569)'} % 209569 is manually annotated
model.grRules{1174} = '(187779)'
model.grRules{1177} = '(274265)' % Manually annotated, ACACT_m
model.grRules{1311} = '(202663)' % Manually annotated, BUT3ZI_m
% model.grRules{1312} = No genes, partial
model.grRules{1313} = '(185524)' % Manually annotated, LIPPT_m
% model.grRules{1351, 1459, 1473, 1484, 1511, 1514, 1516, 1519, 1577, 1579, 1585, 1604, 1643, 1745, 1830} = No genes in P. tricornutum
% model.grRules{1352} = No genes targeted to mitochondria and partial
% sequence
model.grRules{1353} = '(204937)' % Manually annotated, ORPRT_c
model.grRules([1354, 1355:1359], 1) = {'(262343) or (178727)'} 
model.grRules([1360:1364], 1) = {'(267891) or (269435)'} 
model.grRules([1365, 1366], 1) = {'(191760)'}
model.grRules([1367:1378], 1) = {'(211690)'}
% model.grRules{1379 to 1382} = No genes targeted to cytosol
model.grRules{1383} = '(267406)' 
model.grRules{1384} = '(267087) or (268031)' 
model.grRules{1385} = '(170468)' 
model.grRules{1386} = '(170468)'
model.grRules{1387} = '(185726)'
model.grRules{1388} = '(189350) or (215317)' % Based on old network
model.grRules{1389} = '(211024)' % Manually annotated, TMDPP_c
model.grRules{1390} = '(159735)'
% model.grRules{1391, 1392} = No genes, partial
model.grRules{1393} = '(234343)' % Based on the old network
model.grRules{1394} = '(209097)'
model.grRules{1395} = '(167203)' % Based on the old network
model.grRules{1396} = '(167203)'
model.grRules{1397} = '(235828)'
model.grRules{1398} = '(235828)'
model.grRules{1399} = '(229076)'
model.grRules{1400} = '(229076)'
% model.grRules{1401 to 1405} = No genes, partial
model.grRules([1406, 1407], 1) = {'(184309)'} % Manually annotated, HPRT_c
model.grRules([1408:1413], 1) = {'(168304) or (206304)'}
model.grRules([1414:1420], 1) = {'(189350)'}
model.grRules{1421} = '(276698)'  % Old network and blast
model.grRules([1422, 1423], 1) = {'(238041)'}
model.grRules([1424, 1425], 1) = {'(209625)'}  % Manually annotated, HXANDH_c
model.grRules{1426} = '(259696)'
model.grRules{1427} = '(174022)'
% model.grRules{1428} = No significant similarities
% model.grRules{1429} = No genes, partial
model.grRules{1430} = '(205029)'
model.grRules{2155} = '(214292)'
model.grRules([1431:1434], 1) = {'(267891) or (269435)'}
model.grRules([1435:1438], 1) = {'(191760)'}
model.grRules{1439} = '(201468)'
model.grRules{1440} = '(186413)'
model.grRules{1441} = '(222092)' % Based on the old network
model.grRules([1442, 1443], 1) = {'(248379)'}
model.grRules{1444} = '(256661)'
model.grRules{1445} = '(226362) or (178240) or (191087)' 
% model.grRules{1446} = No genes, partial
model.grRules{1447} = '(277202)' % Manually annotated, ME2_m
model.grRules{1448} = '(168262)' % Manually annotated, ICDH_m
model.grRules{1449} = '(274045)'
% model.grRules([1450, 1958], 1) = No genes, partial
model.grRules{2079} = '(247180) or (238236)'
model.grRules{1451} = '(204937)' % Manually annotated, OMPDC_c
% model.grRules{1460} = No significant homologous
model.grRules{1461} = '(174417) or (233742) or (197117) or (274777) or (227544)' % ADDITION OF 3 MANUALLY ANNOTATED ISOZYME IN F. CYLINDRUS (there was 2 proteins in P. tricornutum). In total, all the 5 Manually annotated
model.grRules{1462} = '(191456)'
model.grRules{1463} = '(241346)'
model.grRules{1464} = '(181253) or (224007)'
model.grRules{1465} = '(274844)'  % Manually annotated, ADMDC_c
model.grRules{1466} = '(225230) or (172086)'
model.grRules{1467} = '(235576) or (240044)'
model.grRules{1468} = '(235576) or (240044)'
model.grRules{1469} = '(264866)'
model.grRules{1470} = '(264866)'
model.grRules{1471} = '(178356)'
model.grRules{1472} = '(178356)'
model.grRules{1474} = '(232580 and 269253)'
model.grRules{1475} = '(264525)'
% model.grRules{1476, 1478} = No chloroplastic genes
model.grRules{1477} = '(268504)'
model.grRules{1479} = '(268504)'
model.grRules([1480:1483], 1) = {'(271202)'}
% model.grRules([1485:1486], 1) = No significant similarities
model.grRules{1487} = '(212372)'
model.grRules{1489} = '(186746)'  % Manually annotated, GF6PTA_c
model.grRules{1490} = '(250461)'
model.grRules{1491} = '(246432)'
% model.grRules{1492, 1493} = No significant similarities.
model.grRules{1494} = '(185767)'
% model.grRules{1586} = No genes, partial
% model.grRules{1495} = No significant similarities.
model.grRules{1496} = '(187779)'
model.grRules{1497} = '(187779)'
model.grRules{1498} = '(202663)'
model.grRules{1499} = '(202663)'
model.grRules{1500} = '(223467)'
model.grRules{1501} = '(264261) or (229426)'
model.grRules{1502} = '(229426) or (270026)'  % Manually annotated, 
model.grRules{1503} = '(167310 and 136622)'
model.grRules{1504} = '(271197)'
model.grRules{1505} = '(273786)'
model.grRules{1506} = '(246826)'  % Manually annotated
model.grRules{1508} = '(200742)'
model.grRules{1509} = '(200742)'
% model.grRules{1510} = No mitochondrial genes
model.grRules{1512} = '(178412)'
% model.grRules{1529} = No genes, partial
model.grRules{1513} = '(260099)'
model.grRules{1515} = '(268332)'  % Manually annotated, GDPMEP_c
model.grRules{1517} = '(223160)'  
model.grRules{1518} = '(207230)'  
model.grRules{1520} = '(224930) or (227121)' % Manually annotated, GLACDH_m, addition of one isozyme
% model.grRules{1523} = No cytosolic ascoarbate peroxidase genes
model.grRules{1524} = '(173923) or (152257) or (287495)' % Manually annotated, ASCBPOX_h, addition of two isozymes
model.grRules{1525} = '(287495)'  % Manually annotated
model.grRules{1526} = '(260823) or (196320) or (229440)' % Manually annotated, DHAOX_c
model.grRules{1527} = '(235144)' 
model.grRules{1528} = '(180972)'
model.grRules{1530} = '(233106)'
model.grRules{1533} = '(177880)'
model.grRules{1534} = '(177880)'
model.grRules{1535} = '(269716)'
model.grRules{1536} = '(250461)'
model.grRules{1537} = '(170472)'
model.grRules{1538} = '(217040)'
model.grRules{1539} = '(267900)'
model.grRules{1540} = '(264494) or (291569) or (276598) or (243258) or (246249)' % All 5 are manually annotated CA_h
model.grRules{1541} = '(264424) or (262982) or (260956)' % All are manually annotated CA_m
model.grRules{1547} = '(243875)'  % Manually annotated, NAHCO3t_h
model.grRules{1548} = '(291568)'  % Manually annotated, NAHCO3t_m
model.grRules{1549} = '(259908) or (248163)'
model.grRules{1553} = '(229085)'
model.grRules{1554} = '(229085)'
% model.grRules{1557} = 'No chloroplastic genes
model.grRules{1558} = '(268235)'
model.grRules{1560} = '(268792 and 217881 and 259539 and 268784) or (268792 and 217881 and 259539 and 212206)'
model.grRules{1561} = '(267145 and 268784) or (267145 and 212206)'
model.grRules{1562} = '(154723 and 241864)' % 154723 is partial, but included since this is a potential subunit
model.grRules{1563} = '(273942)'
% model.grRules{1564, 1565} = No genes, partial
model.grRules{1566} = '(183873)' % Manually annotated, G12MT1_c
model.grRules{1567} = '(183873)' % Manually annotated, G12MT2_c
model.grRules{1569} = '(196804)' % Manually annotated, G12MT3_c
model.grRules{1570} = '(196804)' % Manually annotaed, G12MT4_c
model.grRules{1572} = '(235973)'
model.grRules{1573} = '(208572)'
model.grRules{1574} = '(183870) or (187550)'
model.grRules{1575} = '(193959)'
% model.grRules{1576} = No genes, partial
model.grRules{1578} = '(208091) or (169299)' % Manually annotated, DOLASNT_c
model.grRules{1588} = '(168118)' % Manually annotated, MM1_c
model.grRules{1581} = '(189180)'
model.grRules{1583} = '(248439) or (199970)'
model.grRules{1584} = '(170472) or (179924) or (250461)'
model.grRules{1587} = '(188664)'
model.grRules{1588} = '(165775)'
model.grRules{1589} = '(183994)'
model.grRules{1590} = '(205676)' % Manually annotated, MECDPS_h
model.grRules{1596} = '(189402)'
model.grRules{1597} = '(291522)'
model.grRules{1598} = '(242143)'
model.grRules{1599} = '(242143)'
% model.grRules{1601} = No mitochondrial genes
% model.grRules{1603} = No mitochondrial genes
% model.grRules{1605} = No genes, partial
model.grRules{1606} = '(179175 and 224960 and 232303 and 238524)'
model.grRules{1609} = '(291503 and 224960 and 232303 and 238524)'
model.grRules{1611} = '(224960 and 232303 and 238524)'
model.grRules{1612} = '(194850 and 224960 and 232303 and 238524)'
model.grRules{1613} = '(242306 and 224960 and 232303 and 238524)'
model.grRules{1614} = '(291503 and 224960 and 232303 and 238524)'
model.grRules{1615} = '(248177)'
% 1616: No significant similarities.
model.grRules{1617} = '(192704 and 184367 and 276030)'
model.grRules{1618} = '(192704 and 184367 and 276030)'
% 1619: No significant similarities
% 1620 : No genes, partial
model.grRules{1621} = '(226008)'
model.grRules{1622} = '(196266)'
model.grRules{1623} = '(223447)'
model.grRules{1624} = '(223558)'
model.grRules{1625} = '(277329)'
model.grRules{1631} = '(277329)'
% 1634, 1635 : No genes, partial
model.grRules{1626} = '(187842)'
% 1628, 1632, 1633 : No genes, partial
% 1627, 1629 : No significant homologous genes
model.grRules{1630} = '(207407)'
model.grRules{1636} = '(211817)'
model.grRules{1639} = '(225659)'
model.grRules{1846} = '(143441)'
% 1886 : No genes, partial
% 1640 : No genes, partial
model.grRules{1641} = '(205484)'
model.grRules{1642} = '(275015)'
model.grRules{1644} = '(241529)'
model.grRules{1645} = '(241529)'
% Check 1646 : number 500...
% 1647: No chloroplastic genes
% 1648 : No genes, partial (same than 1640)
% 1649, 1679 : No genes, partial
model.grRules{1650} = '(186232)'
model.grRules{1651} = '(186232)'
model.grRules{1652} = '(239507)'
model.grRules{1653} = '(184738)'
% 1655 : No genes, partial (same than 1640)
model.grRules{1656} = '(153813)'
model.grRules{1657} = '(153813)'
model.grRules{1658} = '(196053)'
model.grRules{1659} = '(196053)'
model.grRules{1661} = '(232637)'
model.grRules{1662} = '(224685)'
% 1665 : No mitochondrial genes
model.grRules{2142} = '(238199)'
% 1666, 1667: No mitochondrial genes
model.grRules{1668} = '(269453)'
% 1669: No genes, partial
model.grRules{1671} = '(197036)'
model.grRules{1672} = '(197036)'
% 1673, 1674: No chloroplastic genes
% 1675, 1676, 1677 : No genes, partial
model.grRules{1678} = '(223970)'
model.grRules{1680} = '(259092)'
model.grRules{1681} = '(207169)'
model.grRules{1682} = '(186707)'
model.grRules{1683} = '(267964) or (229079)'
model.grRules{1684} = '(267964) or (229079)'
model.grRules{1685} = '(182520)'
model.grRules{1686} = '(249648)'
model.grRules{1687} = '(264939)'
model.grRules{1688} = '(264939)'
model.grRules{1689} = '(171176)'
model.grRules{1690} = '(206558)'
% model.grRules{1691, 1692} = No genes, partial.
model.grRules{1693} = '(178869) or (225445)'
% 1694, 1695 : No mitochondrial and chloroplastic genes
% 1697, 1698, 1699, 1702, 1706, 1707  : No genes, partial
model.grRules{1704} = '(185157)'
% 1705 : No chloroplastic genes, Selenocysteine-tRNA ligase, chloroplast
model.grRules([1709, 1710, 1711, 1715, 1716, 1717, 1718, 1725, 1726, 1728, 1729, 1731], 1) = {'(264939)'}
model.grRules([1712, 1713, 1719, 1720, 1722], 1) = {'(189189)'}
% 1744, 1749, 1750, 1846, 1871, 1912, 1921, 1987, 1846 : No genes, partial.
model.grRules{1751} = '(192050)'
model.grRules{1747} = '(181873)'
model.grRules{1748} = '(274468)'
model.grRules{1753} = '(291486)'
model.grRules{1754} = '(268562)'
model.grRules{1757} = '(180278)'
model.grRules{1758} = '(196011)'
% 1759, 1760, 1763: No mitochondrial genes
% 1761: No genes, partial
model.grRules{1762} = '(237199)'
model.grRules{1764} = '(180590)'
% 1766: No genes, partial
model.grRules{1767} = '(183660)'
model.grRules{1769} = '(186985)'
model.grRules{1772} = '(244684)'
model.grRules{1773} = '(271325)'
model.grRules{1774} = '(180884) or (241046)'
model.grRules{1775} = '(269134)'  % Manually annotated, PHYEBILIN_h
model.grRules{1779} = '(260963)'
model.grRules{1780} = '(291550)'
model.grRules{1781} = '(291551)'
model.grRules{1782} = '(226929) or (206370) or (225509) or (232063)' % Manually annotated, PLYCOI_h
model.grRules{1783} = '(246826)' % Manually annotated, LYCBC2_h
model.grRules{1784} = '(195646)'
model.grRules{1785} = '(261383)'
model.grRules{1786} = '(232148) or (260743)'  % Manually annotated, ZXANEPOX_h
model.grRules{1787} = '(232148) or (260743)'
model.grRules{1794} = '(232148) or (260743)'
model.grRules{1789} = '(291552) or (212709) or (267113)' % Manually annotated
model.grRules{1790} = '(291552) or (212709) or (267113)'
model.grRules{1793} = '(291552) or (212709) or (267113)'
model.grRules{1797} = '(228459)'
model.grRules{1798} = '(166658)'
% 1799 : No chloroplastic genes
model.grRules{1800} = '(229494)'
% 1801 : No chloroplastic genes
% model.grRules{1802} = No mitochondrial genes
model.grRules{1803} = '(170146)'
model.grRules{1804} = '(205100)'
model.grRules{1805} = '(192916)'
model.grRules([1806, 1807:1811], 1) = {'(260128)'} % Manually annotated, GDR_h
model.grRules{1812} = '(182113) or (185750)' % Manually annotated, GTHP_c
model.grRules([1813, 2147], 1) = {'(275587 and 271949) or (208164 and 271949) or (262343 and 271949)'}
model.grRules{1840} = '(169676)'
model.grRules{1849} = '(232262)' % Manually annotated FK_c
model.grRules{1850} = '(232262)'
model.grRules{1851} = '(206067)' % Manually annotated, GMAND_c
model.grRules{1852} = '(271230)' % Manually annotated, GFUCS_c
model.grRules{1853} = '(263182)' % Manually annotated, same than PFK (rxn 385)
model.grRules{1854} = '(259733) or (254523)' % Manually annotated, RBKRPE_h. Addition of one more isozyme.
model.grRules{1855} = '(183949)'
model.grRules{1856} = '(261602)'
model.grRules{1857} = '(211962)'
model.grRules{1858} = '(264529)'  % Manually annotated, GALK_c
model.grRules{1859} = '(196785)' % Manually annotated, UDPGLDC_c
model.grRules{1860} = '(184505)'
model.grRules{1861} = '(191008) or (170452)' % Manually annotated, UDPRHMS_c
model.grRules{1862} = '(169332)'  % Manually annotated, CPS_m 
model.grRules{1863} = '(264682)'  % Manually annotated, CPS_c
model.grRules{1865} = '(214290) or (212021) or (189282)'
model.grRules{1868} = '(184922) or (238253) or (215191)'
model.grRules{1869} = '(264161)'
model.grRules{1870} = '(264161)'
model.grRules{1872} = '(178244)'
model.grRules{1879} = '(247054)' % Manually annotated,ATPt_h
model.grRules{1880} = '(168527)'
model.grRules{1892} = '(213003)'
model.grRules{1884} = '(277748)' % Manually annotated, MSBENZMT_h
model.grRules{1885} = '(173266)'
model.grRules{1893} = '(193801)'
model.grRules{1898} = '(193801)'
% 1894, 1896 : No genes, partial
% 1899: No chloroplastic genes
model.grRules{1906} = '(212621)'
model.grRules{1907} = '(212621)'
model.grRules{1908} = '(276664)'
model.grRules{1909} = '(276664)'
model.grRules{1911} = '(277202)' % Manually annotated, ME_x
model.grRules{1913} = '(256661)'
model.grRules{1914} = '(268163 and 267319)'
model.grRules{1941} = '(197670)'
model.grRules{1954} = '(253243)'
model.grRules{1955} = '(190546) or (188784)' % Manually annotated,FTR_h
model.grRules{1956} = '(232486)' % Manually annotated, NTRC_h
model.grRules{1959} = '(247054)'
model.grRules{1960} = '(247054)' % Manually annotated, ADPHt_h
model.grRules{1962} = '(224245)'
model.grRules{1971} = '(207890)'
model.grRules{1985} = '(191403)'
model.grRules{1986} = '(213003)'
model.grRules{1990} = '(196011)'
model.grRules{1991} = '(196011)'
model.grRules{1992} = '(196011)'
model.grRules{2001} = '(209569 and 202663 and (270026 or 229426) and 268162)'
model.grRules{2002} = '(195167)'
model.grRules{2003} = '(223740)'
model.grRules{2004} = '(209086)' % Manually annotated, UAGDP_c
model.grRules{2031} = '(209641)' % Manually annotated, GLUCYS_h
% 2039: No significant similarities
% 2041 : No genes, partial
model.grRules{2044} = '(177727) or (181131)'  % Based on old network
model.grRules{2058} = '(168527)'
model.grRules{2059} = '(208164)'
model.grRules{2064} = '(246398)'
model.grRules{2060} = '(226726)'
model.grRules{2061} = '(208572)'
model.grRules{2062} = '(208572)'
% 1579, 2063: No genes in P. tricornutum
model.grRules{2065} = '(267437)'
model.grRules{2066} = '(205516)'
model.grRules{2067} = '(213110)'
model.grRules{2068} = '(226798)'
model.grRules{2069} = '(206154)'
model.grRules{2070} = '(188671)'
% 2071: No genes, partial
% 2072: No significant homologous
model.grRules{2073} = '(261851)'
model.grRules{2074} = '(261851)'
model.grRules{2075} = '(264535) or (273871)' % Both manually annotated, ATPS_c
model.grRules{2076} = '(170931)' % Manually annotated, XYLE_c
model.grRules{2078} = '(236827)'
model.grRules{2103} = '(225753)'
model.grRules{2104} = '(225753)'
model.grRules{2123} = '(209197 and 182880)'
model.grRules{2124} = '(191160)'
model.grRules{2125} = '(264103) or (246693) or (244179) or (291556)' % Addition of 2 other isozymes. All manually annotated
model.grRules{2141} = '(261303)'  % Manually annotated, B13GS
model.grRules{2146} = '(239458)'  % Manually annotated, SOD_m
model.grRules{2155} = '(214292) or (224737)'
model.grRules{1457} = '(259813)'
model.grRules{1458} = '(275677)'
model.grRules{1488} = '(208511)'
model.grRules{1783} = '(246826)' % Manually annotated, LYCBC2_h

% For 22 protein IDs out of 1025 (2.2%) in the P. tricornutum model , we found no homologous proteins.
% For 46 protein IDs (4.5%), partial sequences were found. 
% For 49 proteins (4.8%), predicted subcellular location in F. cylindrus was not the same than in P. tricornutum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Saving a .xml file . Here the sbml file is valid!
% Try keeping a part of the objective function, but set integer coefficients
% for metCharges, metFormulas and coefficients
models = model
MetCharges_Rm = [1658,1659,1660,1662,1663,1668,1669]
MetCharges_Rm = MetCharges_Rm'
i = 0
for i = 1:length(MetCharges_Rm)
    models.metCharges(MetCharges_Rm(i)) = 0
end
models.metCharges([1658,1659,1660,1662,1663,1668,1669])

MetFormulas_Rm = [1658,1659,1660,1661,1662,1663,1668,1669]
MetFormulas_Rm = MetFormulas_Rm'
i = 0
for i = 1:length(MetFormulas_Rm)
    models.metFormulas(MetFormulas_Rm(i)) = {''}
end
models.metFormulas([1658,1659,1660,1661,1662,1663,1668,1669])

Rxn_Bio = {'biomass_pro_c', 'biomass_DNA_c', 'biomass_RNA_c', 'biomass_pigm_h',...
    'biomass_mem_lipids_c', 'biomass_TAG_c', 'biomass_carb_c', 'DM_biomass_c', 'bof_c'}
ID = findRxnIDs(models, Rxn_Bio)
printRxnFormula(models, Rxn_Bio)
checkObjective(models)
ID = findRxnIDs(models, (checkObjective(models)))
Rxn_Bio = {'biomass_pro_c', 'biomass_DNA_c', 'biomass_RNA_c', 'biomass_pigm_h',...
    'biomass_mem_lipids_c', 'biomass_TAG_c', 'biomass_carb_c', 'bof_c'}
ID = findRxnIDs(models, Rxn_Bio)
% Remove biomass equations (except 'DM_biomass_c', the obj function) with non-integer
% coefficients
models = removeRxns(models, Rxn_Bio)

% Writing the model in .xml
writeCbModel(models, 'fileName', 'UpdatedFcylindrus2.xml')

% Save a .xls file (this can saves all the model)
writeCbModel(models, 'format','xls', 'fileName', 'UpdatedFcylindrus.xls')

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updating the model with the suggested modifications of Broddrick et al
% 2019

% Delete 17 reactions
length(model.rxns)
length(model.mets)
Rxn_set_to_remove = {'PGTPt_m'; 'PGTPti_h'; 'GLUDC_m'; 'DHAD_4_m'; '4MOPt_h';...
    '4MOPt_m'; '3MOPt_h'; '3MOPt_m'; 'Ru5PPi_th'; 'Xu5PPi_th'; '6PGTt_h';...
    'AKGMALt_h'; 'MALOAAt_h'; 'MALICITt_m'; 'AKGICITt_m'; 'OAAICITt_m'; 'ME_x'}
model = removeRxns(model, Rxn_set_to_remove)
length(model.rxns)
length(model.mets)

% Move the reaction 'EDD_m' to the plastid , i.e., EDD_h
% But, 2ddg6p_h does not exist
% model1 = addReaction(model1, 'EDD_h', 'metaboliteList', {'6pgc_h',...
%    'h2o_h', '2ddg6p_h'},...
%    'stoichCoeffList', [-1; 1; 1],...  
%    'reversible', true);
% I decided not to add this reaction since it creates a dead-end in the chloroplast
% (2ddg6p_h)

% Move 'ARG_c' to mito 'ARG_m'
% arg__L_m does not exist, must add a transporter from cytosol to mito
ID_ARG_c = findRxnIDs(model, 'ARG_c')
model.rxns(ID_ARG_c) = {'ARG_m'}
model = addReaction(model, 'ARG_m', 'metaboliteList', {'h2o_m',...
    'arg__L_m', 'orn_m', 'urea_m'},...
    'stoichCoeffList', [-1; -1; 1; 1],...  
    'reversible', false);
model.grRules{findRxnIDs(model, 'ARG_m')} = '(182884)'
model.metCharges(findMetIDs(model, 'arg__L_m')) = 1
model.metFormulas(findMetIDs(model, 'arg__L_m')) = {'C6H15N4O2'}

% Add 'ARGt_m'
model = addReaction(model, 'ARGt_m', 'metaboliteList', {'arg__L_c',...
    'arg__L_m'},...
    'stoichCoeffList', [-1; 1],...  
    'reversible', false);
% Gene of ARGt_m : No genes

% Correct the reaction name error 'GLYTA_c' by 'GLYTA_m'
ID_GLYTA = findRxnIDs(model, 'GLYTA_c')
model.rxns(ID_GLYTA) = {'GLYTA_m'}

% Add 4 proteins (SPT_m, HYPRRx_m, GLYCK2_m, ASPGLU2_m)
% 'hpyr_m' does not exist, should be updated
length(model.rxns)
length(model.mets)
model = addReaction(model, 'SPT_m', 'metaboliteList', {'pyr_m',...
    'ser__L_m', 'ala__L_m', 'hpyr_m'},...
    'stoichCoeffList', [-1; -1; 1; 1],...  
    'reversible', true);
% pyr_m + ser__L_m <-> ala__L_m + hpyr_m
model.metCharges(findMetIDs(model, 'hpyr_m')) = -1
model.metFormulas(findMetIDs(model, 'hpyr_m')) = {'C3H3O4'}
% Gene of SPT_m : partial sequences

%{
 test = contains(model1.rxns, 'ASPGLU')
 find(test == 1)
 model1.rxns(ans)
 
 test = contains(model1.mets, 'hpyr')
 find(test == 1)
 model1.mets(ans)
 %}

% 'glyc__R_m' does not exist, should be updated
length(model.rxns)
length(model.mets)
model = addReaction(model, 'HYPRRx_m', 'metaboliteList', {'h_m',...
    'hpyr_m', 'nadh_m', 'glyc__R_m', 'nad_m'},...
    'stoichCoeffList', [-1; -1; -1; 1; 1],...  
    'reversible', true);
% 'h_m + hpyr_m + nadh_m <=> glyc__R_m + nad_m'
model.metCharges(findMetIDs(model, 'glyc__R_m')) = -1
model.metFormulas(findMetIDs(model, 'glyc__R_m')) = {'C3H5O4'}
model.grRules{findRxnIDs(model, 'HYPRRx_m')} = '(267482)'

model = addReaction(model, 'GLYCK2_m', 'metaboliteList', {'atp_m',...
    'glyc__R_m', '2pg_m', 'adp_m', 'h_m'},...
    'stoichCoeffList', [-1; -1; 1; 1; 1],...  
    'reversible', false);
% 'atp_m + glyc__R_m --> 2pg_m + adp_m + h_m'
% Gene of GLYCK2_m : partial sequences

model = addReaction(model, 'ASPGLU2_m', 'metaboliteList', {'asp__L_m',...
    'glu__L_c', 'asp__L_c', 'glu__L_m'},...
    'stoichCoeffList', [-1; -1; 1; 1],...  
    'reversible', true);
% 'asp__L_m + glu__L_c <=> asp__L_c + glu__L_m'
model.grRules{findRxnIDs(model, 'ASPGLU2_m')} = '(189026)'

% Modify Rxn #380 PSI_u and #288 NADHOR_m (complex 1) so that energetic coupling of Bailleul et al is
% included
% 1 h_c and 1 nad_m production by complex 1 (and hence ATP_m) is coupled to 1/6 PSI reaction and hence
%  NADPH generation. So, this should increase ATP:NADPH ratio as CEF_h....
% The more we increase this coupling to 1/5, 1/4... the more the growth
% rate decrease and the more ATP:NADPH ratio increase.
rxns_length1 = length(model.rxns)
mets_length1 = length(model.mets) 
model = addReaction(model, 'PSI_u', 'metaboliteList', {'fdxox_h', 'photon_h', 'pccu1p_u', 'fdxrd_h', 'pccu2p_u', 'mito_coupling_const_h'},...
    'stoichCoeffList', [-2;-2;-2;2;2;-0.15],...  % - 0.15
    'reversible', false);
model = addReaction(model, 'NADHOR_m', 'metaboliteList', {'h_m', 'nadh_m', 'q9_m', 'h_c', 'nad_m', 'q9h2_m', 'mito_coupling_const_h'},...
    'stoichCoeffList', [-5;-1;-1;4;1;1;1],...
    'reversible', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updating the .rules and .genes field using the new .grRule field
%
% This requires updating the model.genes field first.
% Extracting genes from 'grRules' and creating the model.genes field

model.grRules
expression = '\d*\d*'
matchNumber = regexp(model.grRules,expression,'match')
% Unest the cell array
flatCellArray = [matchNumber{:}]
% Transpose the row cell array to a column cell array
flatCellArray = flatCellArray'
% Fill the .genes field with the list of genes
model.genes = flatCellArray
% Removing the same genes occuring more than one time
% By default, the unique function uses the argument 'sorted' and replace (if needed) in
% ascending order the number in the vector. 
model.genes = unique(model.genes)
% generate the .rules field with the .grRules field and the .gene field
model = generateRules(model)
model = buildRxnGeneMat(model) % A matrix with n reactions by m genes

% save modelMichel.mat model
%{
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model, 'FBA');
i =0
for i = 1:10
    disp(i)
    [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model, 'FBA', model.genes(i));
end
%}

%{
% Find the non-empty cells in model.grRules
index = find(~cellfun(@isempty,model.grRules))   
% Fixing problems with genes
model.genes() % This gives a list of old genes
Gene_length = length(index) % Number of non-empty grRules
geneCellarray = cell(Gene_length,1)
% convert a cell array (double) to character
geneCellArrayChar = cellfun(@num2str, geneCellarray, 'UniformOutput', false)
% change model.genes to character
model.genes = geneCellArrayChar

% Creating a cell array of gene numbers using only the model.grRules with
% values
i = 0
for i = 1:length(index)
    model.genes{i} = model.grRules{index(i)} 
end

% Creating a cell array of gene numbers with one gene number per cell
model.genes
expression = '\d*\d*'
matchNumber = regexp(model.genes,expression,'match')
% Unest the cell array
flatCellArray = [matchNumber{:}]
% Transpose the row cell array to a column cell array
flatCellArray = flatCellArray'
% Fill the .genes field with the list of genes
model.genes = flatCellArray
% Removing the same genes occuring more than one time
% By default, the unique function uses the argument 'sorted' and replace (if needed) in
% ascending order the number in the vector. 
model.genes = unique(model.genes)
% generate the .rules field with the .grRules field and the .gene field
model = generateRules(model)

% Do not delete the rxnGeneMat field, which is necessary for the command
% [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(model, 'FBA')
model = buildRxnGeneMat(model)
model = generateRules(model)
verifyModel(model, 'simpleCheck', true)
%}

% Display the number of genes or protein IDs in the model
NumberGenes = length(model.genes);
fprintf('Number of genes or proteins IDs : %4.0f\n',NumberGenes)

% Display the number of reactions
NumRxn = length(model.rxns);
fprintf('The total number of reactions is : %4.0f\n',NumRxn)

% Display the number of Gene-product relationships in the grRules field
% (removing blank entries)
R = length(rmmissing(model.grRules));
fprintf('The total number of GPR is : %4.0f\n',R)

% Printing the number of missing GPR
fprintf('There are %4.0f missing GPR\n', (NumRxn-R))

% Number of reactions with isozymes ('OR' in .grRules)
Index_or = find(contains(model.grRules,'or'));
fprintf('There are %4.0f reactions with isozymes \n', length(Index_or))

% Number of reactions with subunits ('and' in .grRules)
Index_and = find(contains(model.grRules,'and'));
fprintf('There are %4.0f reactions with subunits \n', length(Index_and))

% Number of reactions with only one gene
Numb_One_Gene = length(model.grRules) - length(Index_or) - length(Index_and) - (NumRxn-R);
fprintf('There are %4.0f reactions with only one gene \n', Numb_One_Gene)

% Number of unique GPR
fprintf('There are %4.0f unique GPR \n', length(unique(model.grRules)))

% Number of repeated GPR
fprintf('There are %4.0f repeated GPR \n', R - length(unique(model.grRules)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check and set constraints
% copy model to model1
model1 = model
printConstraints(model1, -1000, 1000)
[selExc, selUpt] = findExcRxns(model1) %, 0, 0)
uptakes = model1.rxns(selUpt) %  Only uptake reactions
exchanges = model1.rxns(selExc) % Uptake, loss, DM, sink reactions

% Fix boundaries of all exchange reactions excepted 'EX_photon_e'
% Assume the N source is NO3 as in laboratory cultures using the f/2 medium.
% hco3_e is the carbon source... NOT co2_e !
EXrxn = {'EX_co2_e','EX_hco3_e','EX_no2_e','EX_no3_e','EX_nh4_e','EX_biotin_e','EX_fe2_e','EX_h_e','EX_h2o_e','EX_o2_e','EX_pi_e','EX_na1_e','EX_so4_e','EX_hso3_e','EX_mg2_e','EX_glyclt_e','EX_selt_e','EX_glc_e','EX_cl_e','EX_thm_e','EX_h2_e','EX_fol_e','EX_co_e','EX_cyan_e','EX_cynt_e','EX_tcynt_e','EX_lac_d_e','EX_etoh_e'}
EXub = [1000,        1000,       0,          0,          1000,     0,            1000,      1000,    1000,     1000,     1000,     0,        0,          0,          0,         0,            0,          0,          0,        0,        0,        0,0,0,0,0,0,0]
EXlb = [0,          -1000,       0,         -1000,       0,       -1000,         -1000,     -1000,   -1000,    -1000,    -1000,   -1000,    -1000,       0,        -1000,       0,            0,          0,        -1000,     -1000,     0,        0,0,0,0,0,0,0]
for i = 1:length(EXrxn);
    index = findRxnIDs(model1,EXrxn{i});
    model1.ub(index) = EXub(i);
    model1.lb(index) = EXlb(i);
end

% No constraints on CEF, it can be calculated in the simulation
model1 = changeRxnBounds(model1,'CEF_h', 0,'l');
model1 = changeRxnBounds(model1,'CEF_h', 1000,'u');

% ATPt_h  was set with boundaries = 0 (lower and upper!) I changed it
model1 = changeRxnBounds(model1,'ATPt_h',0,'l');
model1 = changeRxnBounds(model1,'ATPt_h',1000,'u');

% We can prevent futile cycle from carrying any flux although those
% reactions are probably inactive with our two-step FBA
rxns = {'DIADINXDE_h','VIOXANDE_h','ANTHXANDE_h'};
for r = 1:length(rxns)
    model1 = changeRxnBounds(model1,rxns{r},0,'b');
end

% Set reaction fluxes of ATPM_c, ATPM_m, ATPM_h = 0...
% NGAM is set with AOX_m . Assume xx mmol / g DW / h (see below)
model1 = changeRxnBounds(model1, 'ATPM_c', 0, 'b')
model1 = changeRxnBounds(model1, 'ATPM_m', 0, 'b')
model1 = changeRxnBounds(model1, 'ATPM_h', 0, 'b')

% Set boundaries of 'PDYXPT and PYDXDH ' = 0, otherwise an unrealistic flux
% loop is posible
model1 = changeRxnBounds(model1, 'PDYXPT_c', 0, 'b')
model1 = changeRxnBounds(model1, 'PYDXDH_c', 0, 'b')

% Set PTOX and MEHLER reactions = 0 as in P. triconrutum (Bailleul et al 2015)
% Respiration and Mehler reaction are thought to be low in polar diatoms
% (Goldman et al )
model1 = changeRxnBounds(model1, 'PTOX_h', 0, 'b')
model1 = changeRxnBounds(model1, 'MEHLER_h', 0, 'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%0.5 g C / g DW = 0.0417 mol C / g DW ou 41.7 mmol C / gDW 
%41.7 * (0.4/24) = 0.69 mmol C / g DW / h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updating the biomass objective function

% Pigments measured in F. cylindrus under continuous light (30 umol m-2 s-1) (Guérin et al, in prep)
Ccell = 7.5E-12   % g C / cell
DWcell = Ccell * 2. % 1.5E-11 % Ccell * 2.  % Assume a g DW / g C = 2, Ccell * 2. = 1.5e-11
chla_P = 0.14          % pg Chla / cell
chlc_P = 0.03           % pg Chlc / cell
Fucoxan_P = 0.06        % pg Fucoxanthin / cell
B_carot_P = 0.004       % pg B-carotene / cell
DD_P = 9.0/1000        % pg DD / cell   (Alderkamp et al., 2011)

% Pigments in percentage of mass per g DW (% g / gDW)
chla = chla_P * 1E-12 / DWcell * 100.; 
Fucoxan = Fucoxan_P * 1E-12 / DWcell * 100.;
B_carot = B_carot_P * 1E-12 / DWcell * 100.;
DD = DD_P * 1E-12 / DWcell * 100.;
chlc = chlc_P * 1E-12 / DWcell * 100.;
chlc_rat = 2.2;    % Ratio chlorophyll c2 : chlorophyl c1 in P. tricornutum (Owens & Wold (1986) Plant Physiol 80: 732?738)		
chlc1 = 1/chlc_rat * chlc_P;
chlc2 = (1-(1/chlc_rat)) * chlc_P;

% Molar mass of components (g/mol)
molChla =	892.4782;
molFucoxan =	658.9040;
molB_carot =	536.8704;
molDD = 582.8528;
molChlc1 =	607.9166;
molChlc2 =	609.9324;

% Pigment composition (mol g DW-1)
chla = chla / molChla;
Fucoxan = Fucoxan / molFucoxan;
B_carot = B_carot / molB_carot;
DD = DD / molDD;
chlc1 = chlc1 / molChlc1;
chlc2 = chlc2 / molChlc2;

% Non-normalized biomass weight (g / g DW)
NonNormBiom = (chla * molChla) + (Fucoxan * molFucoxan) + (B_carot * molB_carot) + (DD * molDD) + (chlc1 * molChlc1) + (chlc2 * molChlc2);
Correction = NonNormBiom; 

% Stoichiometric coefficients (mmol / g DW)
chla = chla / Correction * 1000.;
Fucoxan = Fucoxan / Correction * 1000.;
B_carot = B_carot / Correction * 1000.;
DD = DD / Correction * 1000.;
chlc1 = chlc1 / Correction * 1000.;
chlc2 = chlc2 / Correction * 1000.;

printRxnFormula(model1, 'biomass_pigm_h');
model1 = addReaction(model1, 'biomass_pigm_h', 'metaboliteList', {'cholphya_h',...
    'cholphyc1_h', 'cholphyc2_h', 'caro_h', 'fxanth_h',...
    'diadinx_h', 'biomass_pigm_h'},...
    'stoichCoeffList', [-chla; -chlc1; -chlc2; -B_carot; -Fucoxan;...
    -DD; 1],...
    'reversible', false);
%rxns_length1 = length(model1.rxns);
%mets_length1 = length(model1.mets); 

% Update biomass_DNA_c for the total F. cylindrus genome size and its low
% G:C content (Mock et al 2017)
GC_P = 0.398    % Proportion of total nucleotides as GC (Mock et al 2017)
AT = 1-GC_P;
GenomeSize_P = 61.1E+6   % Total genome size (number of bases)

% Calculating total number of nucleotides per cell (dAMP, dTMP, dGMP, dCMP)
dNucl = [GenomeSize_P * AT / 2, GenomeSize_P * AT / 2, GenomeSize_P * GC_P / 2, GenomeSize_P * GC_P / 2];

% Molar mass of nucleotides (mol/g) (moldAMP, moldTMP, moldGMP, moldCMP)
molMassNucl = [329.2055, 320.1921, 345.2049, 305.1808];

% Nucleotides (g / cell)
Avog = 6.02E+23;     % Avogadro constant in atom / mol
NuclMass = dNucl .* molMassNucl / Avog;

% Relative nucleotide and pyrophosphate abundance (number or mol/gDW)
Rel_dNucl = dNucl ./ sum(dNucl);
Rel_Pyr = 1;

% triphosphate nucleotide (dNTP) and pyrophosphate molar mass (g/mol)
% (dATP, dTTP, dGTP, dCTP, pyrophosphate or ppi)
MolMassdNTP = [487.1495, 478.1361, 503.1489, 463.1248];
MolMassPyr = 174.9513;

% Non-normalized biomass weight (g / gDW = gDW / gDW)
NonNormBiom = (sum(Rel_dNucl .* MolMassdNTP)) - (Rel_Pyr * MolMassPyr);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients (mmol / gDW) of dATP, dTTP, dGTP, dCTP (S_dNTP) and Pyr
% (S_Pyr)
dNTP = Rel_dNucl ./ Correction;
Pyr = Rel_Pyr / Correction;

printRxnFormula(model1, 'biomass_DNA_c');
model1 = addReaction(model1, 'biomass_DNA_c', 'metaboliteList', {'dgtp_c',...
    'dttp_c', 'dctp_c', 'datp_c', 'ppi_c', 'biomass_dna_c'},...
    'stoichCoeffList', [-dNTP(3); -dNTP(2); -dNTP(4); -dNTP(1); Pyr; 1],...
    'reversible', false);
%rxns_length1 = length(model1.rxns);
%mets_length1 = length(model1.mets);

% Update 'biomass_RNA_c' assuming a RNA:DNA ratio of 8 (Lourenco et al (1998) J Phycol 34: 798?811)		
% RNA Nucleotides (g / cell) (dAMP, dTMP, dGMP, dCMP)
RNA_DNA_r_P = 8.
Mass_Nucl = NuclMass * RNA_DNA_r_P;

% Relative abundance is the same than for DNA
Rel_Nucl = Rel_dNucl;
Rel_Pyr;

% Molar mass of RNA nucleotides (g/mol) ATP, UTP, GTP, CTP, Pyr
MolMassNTP = [503.1489, 480.1090, 519.1483, 479.1242];
MolMassPyr = 174.9513;

% Non-normalized biomass weight
NonNormBiom = (sum(Rel_Nucl .* MolMassNTP)) - (Rel_Pyr * MolMassPyr);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients (ATP, UTP, GTP, CTP, Pyr)
NTP = Rel_Nucl ./ Correction;
Pyr = Rel_Pyr / Correction;

printRxnFormula(model1, 'biomass_RNA_c');
model1 = addReaction(model1, 'biomass_RNA_c', 'metaboliteList', {'atp_c',...
    'gtp_c', 'utp_c', 'ctp_c', 'ppi_c', 'biomass_rna_c'},...
    'stoichCoeffList', [-NTP(1); -NTP(3); -NTP(2); -NTP(4); Pyr; 1],...
    'reversible', false);
%rxns_length1 = length(model1.rxns);
%mets_length1 = length(model1.mets); 

% Update biomass_mem_lipids_c equation
% Molar mass of lipids in g/mol ()
molMassLip = [768, 722, 724, 726, 936, 765, 903, 793, 1123, 764.97,...
    780.97, 824.97, 728.97, 803, 782.97, 806.97, 616, 712, 568]

% Molar ratios of lipids (mol/ gDW) (values for P. tricornutum, Levering et al., 2016)
MolRatLip_P = [0.0555, 0.1244, 0.0622, 0.0504, 0.0450, 0.1043, 0.0318, 0.0144,...
     0.0070, 0.0490, 0.0044, 0.0622, 0.0670, 0.0380, 0.0180, 0.0220,...
     0.0104, 0.0062, 0.0088];

% Mass ratio (g / gDW = gDW / gDW)
MassRatLip = molMassLip .* MolRatLip_P;

% Non-normalized biomass weight (g / gDW = gDW / gDW)
NonNormBiom = sum(MassRatLip);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients of each lipid (mmol / gDW)
S_Lip = MolRatLip_P ./ Correction;

printRxnFormula(model1, 'biomass_mem_lipids_c');
model1 = addReaction(model1, 'biomass_mem_lipids_c', 'metaboliteList', {'mgdg205n3164n1_h',...
    'mgdg1619Z163n4_h', 'mgdg1619Z162n4_h', 'mgdg161_h', 'dgdg205n31619Z_h',...
    'sqdg140160_h', 'sqdg1619Z240_c', 'sqdg160_h', 'asq205n31602o205n3_h',...
    'pg205n31613E_h', 'pc182_9_12_c', 'pc205n3_c', 'pc161_c',...
    'dgta205n3_c', 'pe205n3_c', 'pail1619Z160_c', '12dgr182_9_12_c',...
    '12dgr226n3_c', '12dgr160_h', 'biomass_mem_lipids_c'},...
    'stoichCoeffList', [-S_Lip(1); -S_Lip(2); -S_Lip(3); -S_Lip(4); -S_Lip(5);...
    -S_Lip(6); -S_Lip(7); -S_Lip(8); -S_Lip(9); -S_Lip(10); -S_Lip(11); -S_Lip(12);...
    -S_Lip(13); -S_Lip(14); -S_Lip(15); -S_Lip(16); -S_Lip(17); -S_Lip(18); -S_Lip(19); 1],...
    'reversible', false);
%rxns_length1 = length(model1.rxns);
%mets_length1 = length(model1.mets); 

% Update the 'biomass_pro_c' equation

% Weight amino acid percentage (Brown, 1991). The sum = 100%
% Amino acid list : 'alatrna_c, argtrna_c, asntrna_c, asptrna_c, cystrna_c,
% glntrna_c, glutrna_c, glytrna_c, histrna_c, iletrna_c, leutrna_c,
% lystrna_c, mettrna_c, phetrna_c, protrna_c, sertrna_c, thrtrna_c, trptrna_c,
% tyrtrna_c, valtrna_c
MassPercent_P = [7.2, 6.6, 0, 8.6, 0.38, 0, 11.2, 5.8, 1.6, 4.9, 7.7, 5.6,...
      1.9, 6.6, 6.3, 5.9, 5.4, 1.6, 4.1, 5.6];

% Normalized weight (% g)
NormWeight = MassPercent_P / sum(MassPercent_P); % this makes the total = 1
NormWeight(3) =  MassPercent_P(4) / sum(MassPercent_P) / 2.;
NormWeight(6) = MassPercent_P(7) / sum(MassPercent_P) / 2.;

% Amino acid molar mass (mol/g)
MolMassAmino = [73.0935, 159.2089, 116.1182, 116.0951, 105.1585, 130.1447,...
    130.1216, 59.0670, 139.1548, 115.1730, 115.1730, 131.1955, 133.2115,...
    149.1893, 99.1307, 89.0929, 103.1194, 188.2253, 165.1887, 101.1465];

% Molar Mass of other reactants (g/mol) (Water, atp_c, gtp_c)
MolMassReac = [18.0152, 503.1489, 519.1483];

% Molar Mass of other products (g/mol) (Proton, Phosphate, GDP, ADP,
% Unloaded tRNA)
MolMassProd = [1.0079, 95.9793, 440.1763, 424.1769, 1.0079];

% Molar ratio (% mol/gDW)
MolRatAmino = NormWeight ./ MolMassAmino * 100.;

% Molar ratio of other reactants (mol/ gDW) (Water, atp_c, gtp_c)
SumMolRat = sum(MolRatAmino);
MolRatReac = [3, 1, 2] * SumMolRat;

% Molar ratio other product (mol/ gDW) (Proton, Phosphate, GDP, ADP,
% Unloaded tRNA) where unloaded tRNA accounts for the sum of all molar
% ratio of all amino acids
MolRatProd = [4, 3, 2, 1, 1] * SumMolRat;

% Non-normalized biomass weight (g DW / gDW)
NonNormBiom = sum(NormWeight) * 100. + (sum(MolMassReac .* MolRatReac)) - (sum(MolMassProd .* MolRatProd));
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients of amino acids (mmol / gDW)
S_Amino = - MolRatAmino / Correction;
Sum_S_Amino = - sum(S_Amino);

% Stoichiometric coefficients of other reactants (mmol / gDW)
% Water, ATP, GTP
S_Reac = - [3, 1, 2] * Sum_S_Amino;

% Stoichiometric coefficients of other products (mmol / gDW)
S_Prod = [4, 3, 2, 1] * Sum_S_Amino;

% Stoichiometric coefficients (mmol / gDW) of all 'trna-amino acid species' as products
% 'trnaala_c', 'trnaarg_c'
% But those coefficients are positive...
S_Amino_Prod = - S_Amino;

printRxnFormula(model1, 'biomass_pro_c');
model1 = addReaction(model1, 'biomass_pro_c', 'metaboliteList', {'h2o_c',...
    'atp_c', 'gtp_c', 'alatrna_c', 'argtrna_c', 'asntrna_c',...
    'asptrna_c', 'cystrna_c', 'glntrna_c', 'glutrna_c',...
    'glytrna_c', 'histrna_c', 'iletrna_c', 'leutrna_c', 'lystrna_c',...
    'mettrna_c', 'phetrna_c', 'protrna_c', 'sertrna_c', 'thrtrna_c',...
    'trptrna_c', 'tyrtrna_c', 'valtrna_c',...  % End of reactants
    'h_c', 'pi_c', 'adp_c', 'gdp_c', 'trnaala_c', 'trnaarg_c', 'trnaasn_c',...
    'trnaasp_c', 'trnacys_c', 'trnagln_c', 'trnaglu_c', 'trnagly_c',...
    'trnahis_c', 'trnaile_c', 'trnaleu_c', 'trnalys_c', 'trnamet_c',...
    'trnaphe_c', 'trnapro_c', 'trnaser_c', 'trnathr_c', 'trnatrp_c',...
    'trnatyr_c', 'trnaval_c', 'biomass_pro_c'},...
    'stoichCoeffList', [S_Reac(1); S_Reac(2); S_Reac(3); S_Amino(1); S_Amino(2); S_Amino(3); S_Amino(4);...
    S_Amino(5); S_Amino(6); S_Amino(7); S_Amino(8); S_Amino(9); S_Amino(10); S_Amino(11);...
    S_Amino(12); S_Amino(13); S_Amino(14); S_Amino(15); S_Amino(16); S_Amino(17); S_Amino(18);...
    S_Amino(19); S_Amino(20); S_Prod(1); S_Prod(2); S_Prod(3); S_Prod(4); S_Amino_Prod(1); S_Amino_Prod(2);...
    S_Amino_Prod(3); S_Amino_Prod(4); S_Amino_Prod(5); S_Amino_Prod(6);...
    S_Amino_Prod(7); S_Amino_Prod(8); S_Amino_Prod(9); S_Amino_Prod(10);...
    S_Amino_Prod(11); S_Amino_Prod(12); S_Amino_Prod(13); S_Amino_Prod(14);...
    S_Amino_Prod(15); S_Amino_Prod(16); S_Amino_Prod(17); S_Amino_Prod(18);...
    S_Amino_Prod(19); S_Amino_Prod(20); 1],...
    'reversible', false);
%rxns_length1 = length(model1.rxns);
%mets_length1 = length(model1.mets);


%{
printFluxVector(model1, FBAsolution.x, false, true)
surfNet(model1, 'EX_pi_e', 0, FBAsolution.x)
%}
%{
C_N = 5.72 
% 7.5 or 5.5; Guérin et al measured a C:N ratio of 5.5, 
% Garcia et al measured a C:N ratio of 5.72
N_P = 10.1  % Garcia et al (2018): 10.1 ; N. Schiffrin: N:P= 13
fold_change = 1.1
Ps = 1/12 * (1000 * chla_P * 1E-12 / DWcell) 
N_low =  - Ps / (C_N * fold_change) ; % in mmol / gDW / h % 6.6
N_high = - Ps / (C_N / fold_change);
P_low = - Ps / (C_N * N_P * fold_change);  % in mmol / gDW / h % 106
P_high = - Ps / (C_N*N_P / fold_change); % not 2
%Pho_low = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) * fold_change ; % not 5
%Pho_high = FBAsolution.x(findRxnIDs(model1, 'EX_photon_e')) / fold_change ;

model1 = changeRxnBounds(model1,'EX_no3_e', N_high, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_no3_e', N_low, 'u'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_pi_e', P_high, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_pi_e', P_low, 'u'); % Bd: -1000 / 0
%model1 = changeRxnBounds(model1,'EX_photon_e', Pho_low, 'l'); % Bd: -1000 / 0
%model1 = changeRxnBounds(model1,'EX_photon_e', P_low, 'u'); % Bd: -1000 / 0
FBAsolution = optimizeCbModel(model1, 'max','one')
printFluxVector(model1, FBAsolution.x, false, true)
Interesting results, En forçant le modèle à un ratio C:N de 7.5 et N: P de 13, 
je trouve aucun stockage de NO3- dans la vacuole et un stockage faible de phosphore (environ 1-2% du P intracelllulaire) dans la vacuole. 
Toutefois, quand je fais rouler le modèle avec un ratio de C:N de 5.5 et N: P de 13, 
je trouve un grand (25% du N intracellulaire) stockage de NO3- dans la vacuole , mais pas de stockage de P dans la vacuole !!!
C'est pour une fraction de protéine de 30./100, carb = 20./100,
lipid=10./100
%}


% Add the main objective function with explicit storage_carbon
% (bof_c_accumulation_c)
% Here I replace 'biomass_tag' by 'carbon_storage_c'
% I replaced 'biomass_c' by 'biomass_c_acc_c'

total_prot_P = 40./100 %52./100      %     ratio in g/gDW (selected so that C:N = 6)
total_carb_P = 25./100       %    ratio in g/DW (selected so that C:N = 6)
total_lipid_P = 15./100     %    ratio in g/ gDW (selceted so that C:N = 6)
total_DNA_P = 0.0022        %    ratio in g/gDW calculated for F. cylindrus
total_RNA_P = 0.0176         %    ratio in g/gDW calculated for F. cylindrus
total_pigm_P = 1.56 / 100    %   ratio in g/gDW measured in F. cylindrus
SumTot = total_prot_P + total_carb_P + total_lipid_P + total_DNA_P + total_RNA_P + total_pigm_P

total_prot_Per = total_prot_P * 100.       %     ratio in % g/gDW
total_carb_Per = total_carb_P * 100.      %    ratio in % g/DW
total_lipid_Per = total_lipid_P * 100.     %    ratio in % g/ gDW
total_DNA_Per = total_DNA_P * 100.        %    ratio in % g/gDW calculated for F. cylindrus
total_RNA_Per = total_RNA_P * 100.      %    ratio in % g/gDW calculated for F. cylindrus
total_pigm_Per = total_pigm_P * 100.    %   ratio in % g/gDW measured in F. cylindrus
SumTot = total_prot_Per + total_carb_Per + total_lipid_Per + total_DNA_Per + total_RNA_Per + total_pigm_Per

% corection of main component is necessary because of error in measurements
% and because of unknown components (such as osmolyte) in the cell
% characterization
%prot_cor = (1.-total_pigm_P - total_DNA_P - total_RNA_P) * total_prot_P / (total_prot_P + total_carb_P + total_lipid_P);
%carb_cor = (1. - total_pigm_P - total_DNA_P - total_RNA_P) * total_carb_P / (total_prot_P + total_carb_P + total_lipid_P);
%lip_cor = (1. - total_pigm_P - total_DNA_P - total_RNA_P) * total_lipid_P / (total_prot_P + total_carb_P + total_lipid_P);
prot_cor = ( total_prot_Per + ( 100. - SumTot )/ 3. ) / 100.
carb_cor = ( total_carb_Per + ( 100. - SumTot )/ 3. ) / 100.
lip_cor = ( total_lipid_Per + ( 100. - SumTot )/ 3. ) / 100.
total_DNA = total_DNA_Per / 100.
total_RNA = total_RNA_Per / 100.
total_pigm = total_pigm_Per / 100.
SumTot2 = prot_cor + carb_cor + lip_cor + total_DNA + total_RNA + total_pigm

Prop_glucan = 3*10.      % Proportion of cabohydrate accumulated as glucan
Prop_TAG = 2*10.          % Proportion of lipid accumulated as TAG
mem_lip = lip_cor - (lip_cor * (Prop_TAG/100.))
carb_cor2 = carb_cor - (carb_cor * (Prop_glucan/100.))
TAG = lip_cor * (Prop_TAG/100.);
Glucan = carb_cor * (Prop_glucan/100.);
C_storage = TAG + Glucan;
Mass_bal = prot_cor + carb_cor2 + mem_lip + TAG + Glucan + total_DNA + total_RNA + total_pigm
if Mass_bal < 1.01 & Mass_bal > 0.99;
    disp('Mass balance is ok');
else
    error('Mass of all components is unbalanced');   
end

%rxns_length1 = length(model1.rxns);
%mets_length1 = length(model1.mets); 
model1 = addReaction(model1, 'bof_c_accumulation_c', 'metaboliteList', {'biomass_pro_c',...
    'biomass_pigm_h', 'biomass_mem_lipids_c', 'biomass_dna_c',...
    'biomass_rna_c', 'biomass_carb_c', 'carbon_storage_c', 'biomass_c_acc_c'},...
    'stoichCoeffList', [-prot_cor; -total_pigm; -mem_lip; -total_DNA;...
    -total_RNA; -carb_cor2; -C_storage; 1],...
    'reversible', false);

% Add a reaction 'biomass_c_storage_c'
% Molar ratios (mmol / g DW) (13glucan_c, tag1619Z1619Z160)
% amount accumulate per time step...
%MolRat = [0.117, 0.05]

% Molar mass (g/mol)
Molmass = [419.97, 802];

% Mass ratio (g / gDW)
%MassRat = Molmass .* MolRat / 1000.

% Here we use Glucan and TAG calculated in the mass balance of the
% objective function
MassRat = [Glucan, TAG];

MolRat = MassRat ./ Molmass * 1000; % mmol / gDW

% NonNormalized biomass weight (g / gDW)
NonNormBio = sum(MassRat);
Correction = NonNormBio;

% Stoichiometric coefficients (mmol / g)
S_Stor = MolRat / Correction;

model1 = addReaction(model1, 'biomass_c_storage_c', 'metaboliteList', {'13glucan_c',...
    'tag1619Z1619Z160_c', 'carbon_storage_c'},...
    'stoichCoeffList', [-S_Stor(1); -S_Stor(2); 1],... 
    'reversible', false);

% Update the biomass_carb_c equation

% Molar ratio (mol/ gDW) of the 9 sugar species (Levering et al, 2016)
% Glucose, galactose, mannose, Xylulose, arabinose, fucose, rhamnose,
% glucuronic acid, mannose-sulfate. The total proportion is corrected so
% that all species molar fraction = 1 so that mass balance is conserved.
% sugar_mol = [0.2884	0.1550 0.2475 0.0806 0.0340 0.0436 0.0794 0.0565 0.0201];

% This represents a first approximation, I could have deleted glucose since
% glucan is incorporated explicitly in the model, but this modification
% does not change the modeled growth rate because all species molar fraction is
% reported to 1 so that mass balance is conserved.
% Here we assume negligible free glucose inside the cell.
sugar_mol = [0	0.2206	0.3268	0.1155	0.0447	0.0647	0.1085	0.1005	0.0239];
% Molar mass of sugars (g/mol)
MolMassSug = [564.2851 564.2851 603.3244 534.2592 534.2592 587.3250 548.2857 577.2608 682.3797];

% Mass ratio of sugars (g / g)
MassRatSug = sugar_mol .* MolMassSug;

% Product molar ratio (mol / gDW); 'Proton, UDP, GDP'
Proton_mol = sum(sugar_mol);
UDP_mol = sum(sugar_mol([1,2,4,5,7,8]));
GDP_mol = sum(sugar_mol([3,6,9]));
Prod_mol = [Proton_mol, UDP_mol, GDP_mol];

% Molar mass products (g/mol)
MolMassProd = [1.0079, 401.1370, 440.1763];

% Mass ratio products (g / g)
MassRatSug2 = Prod_mol .* MolMassProd;

% Non-normalized biomass weight (g / gDW)
NonNormBiom = sum(MassRatSug) - sum(MassRatSug2);
Correction = NonNormBiom / 1000.;

% Stoichiometric coefficients (mmol / gDW)
Stoich_Sug = [sugar_mol, Prod_mol];
S_Sug = Stoich_Sug / Correction;

printRxnFormula(model1, 'biomass_carb_c');

model1 = addReaction(model1, 'biomass_carb_c', 'metaboliteList', {'udpg_c',...
    'udpgal_c', 'gdpmann_c', 'udpxyl_c',...
    'udparab_c', 'gdpfuc_c', 'udprmn_c', 'udpglcur_c', 'gdpman2s_c',...
    'h_c', 'udp_c', 'gdp_c', 'biomass_carb_c'},...
    'stoichCoeffList', [-S_Sug(1); -S_Sug(2); -S_Sug(3); -S_Sug(4); -S_Sug(5);...
    -S_Sug(6); -S_Sug(7); -S_Sug(8); -S_Sug(9); S_Sug(10); S_Sug(11);...
    S_Sug(12); 1],...
    'reversible', false);
%rxns_length1 = length(model1.rxns);
%mets_length1 = length(model1.mets);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change Objective function and set constraints

model1 = removeRxns(model1, 'bof_c')
model1 = removeRxns(model1, 'DM_biomass_c')
model1 = removeRxns(model1, 'biomass_TAG_c')

% Add a demand reaction without the command 'addDemand()' otherwise there
% is a bug with the subsystem field...
model1 = addReaction(model1, 'DM_biomass_c_acc_c', 'metaboliteList', {'biomass_c_acc_c'},...
    'stoichCoeffList', [-1],...
    'reversible', false);
model1 = changeObjective(model1, 'DM_biomass_c_acc_c');

Ps = 1/12 * (1000 * chla_P * 1E-12 / DWcell)   % Ps measured by S. Guérin (1 mg C mg chla-1 h-1 / 12 g/mol * mg chla / gDW ) = (mmol C g DW-1 h-1)
% Ps at 30 umol m-2 s-1 at 0 C
% Here I should calculate error propagation on Ps (due to chla, C quotas
% and other parameters...)
%model1 = changeRxnBounds(model1,'EX_no3_e', -Ps/6, 'l'); 
%model1 = changeRxnBounds(model1,'EX_no3_e', -Ps/6, 'u'); 
model1 = changeRxnBounds(model1, 'EX_hco3_e', -Ps, 'b'); % m
model1 = changeRxnBounds(model1, 'EX_co2_e', 0, 'l'); % mmol hco3 g DW-1 h-1
model1 = changeRxnBounds(model1, 'EX_co2_e', 1000, 'u'); 


% Constraints C:N and N:P ratios
C_N = 5.72 
% 7.5 or 5.5; Guérin et al measured a C:N ratio of 5.5, 
% Garcia et al (2018) measured a C:N ratio of 5.72
N_P = 10.1  % Garcia et al (2018): 10.1 ; N. Schiffrin: N:P= 13
fold_change = 1.0
N_low =  - Ps / (C_N * fold_change) ; % in mmol / gDW / h % 6.6
N_high = - Ps / (C_N / fold_change);
P_low = - Ps / (C_N * N_P * fold_change);  % in mmol / gDW / h % 106
P_high = - Ps / (C_N * N_P / fold_change); % not 2


model1 = changeRxnBounds(model1,'EX_no3_e', N_high, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_no3_e', N_low, 'u'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_pi_e', P_high, 'l'); % Bd: -1000 / 0
model1 = changeRxnBounds(model1,'EX_pi_e', P_low, 'u'); % Bd: -1000 / 0

% Set reaction fluxes of ATPM_c, ATPM_m, ATPM_h = ? if needed
model1 = changeRxnBounds(model1, 'ATPM_c', 0.01, 'b')
model1 = changeRxnBounds(model1, 'ATPM_m', 0.01, 'b')
model1 = changeRxnBounds(model1, 'ATPM_h', 0.01, 'b')

model1 = changeRxnBounds(model1, 'EX_photon_e', -1000, 'l'); 
model1 = changeRxnBounds(model1, 'EX_photon_e', 0, 'u');

model1 = changeRxnBounds(model1, 'AOX_m', 0, 'b');

% Perform FBA
FBAsolution = optimizeCbModel(model1, 'max','one')
printFluxVector(model1, FBAsolution.x, false, true)
printFluxVector(model1, FBAsolution.x, true);
uc = FBAsolution.f * 24;

% Printing the predicted growth rate
fprintf('The predicted growth rate is %2.4f d-1 \n', uc)

rxnID_ATPS_m = findRxnIDs(model1, 'ATPS_m');
sol1 = FBAsolution.x(rxnID_ATPS_m);
rxnID_ATPS_h = findRxnIDs(model1, 'ATPS_h');
sol2 = FBAsolution.x(rxnID_ATPS_h);
ratio = sol1/sol2 ; % This is the ratio of mitochondrial ATP produced / chloroplastic ATP produced
fprintf('\nProportion of mitochondrial ATP production: %4.2f \n', sol1 / (sol1+sol2));

% Printing the C:N and N:P molar ratio if no C, N and P loss from the cells
% Be careful, 'EX_' reactions are net fluxes for a given species, but other
% species with C or N might also contribute to uptake/loss.
C_N = FBAsolution.x(findRxnIDs(model1, 'EX_hco3_e')) / FBAsolution.x(findRxnIDs(model1, 'EX_no3_e'));
N_P = FBAsolution.x(findRxnIDs(model1, 'EX_no3_e')) / FBAsolution.x(findRxnIDs(model1, 'EX_pi_e'));
fprintf('C:N molar ratio : %4.2f \n', C_N)
fprintf('N:P molar ratio : %4.2f \n', N_P)

fprintf('\nGlycolate exchange flux : %4.2f \n', FBAsolution.x(findRxnIDs(model1, 'EX_glyclt_e'))) % 0
fprintf('Plastid terminal oxidase flux : %4.2f \n', FBAsolution.x(findRxnIDs(model1, 'PTOX_h'))) % 0
fprintf('Flux through Mehler reaction : %4.2f \n', FBAsolution.x(findRxnIDs(model1, 'MEHLER_h'))) % 0
fprintf('Reverse Rubsico flux or photorespiration : %4.2f \n', FBAsolution.x(findRxnIDs(model1, 'RUBISO_h'))) % photorespiration = 0
fprintf('Production of phosphoglycolate : %4.2f \n', FBAsolution.x(findRxnIDs(model1, 'PGLYCP_h'))) % photorespiration = 0

fprintf('\nWhen using a mitochondrial / chloroplastic coupling constant \n')
fprintf('Flux through AOX : %4.2f \n', FBAsolution.x(findRxnIDs(model1, 'AOX_m'))) % 0 if no coupling constant, or  
fprintf('Cyclic electron flux : %4.2f \n', FBAsolution.x(findRxnIDs(model1, 'CEF_h'))) % 0.2744 if no coupling constant, or 0 with coupling constant

% Calculating proportion of active and inactive reactions
FBAsolution = optimizeCbModel(model1, 'max', 'one');
IndActive = find(FBAsolution.x ~= 0); 
%IndActive = find(abs(FBAsolution.x) > 1E-09);
IndInac = find(FBAsolution.x == 0); 
assert(length(IndActive) + length(IndInac) == length(FBAsolution.x))
fprintf('Number of active reactions : %4.2f \n', length(IndActive) );
fprintf('Proportion of active reactions : %4.2f \n', 100 * length(IndActive) / length(FBAsolution.x));
fprintf('Number of inactive reactions : %4.2f \n', length(IndInac) );
fprintf('Proportion of inactive reactions : %4.2f \n', 100 * length(IndInac) / length(FBAsolution.x));

%{
% Ratio C:N : What ratios should we choose?
N. Schiffrine measured a total C:N ratio with CHN equals to around 7.5 in
cultures of C. gelidus with 100 uM total NO3-, but a C:N ratio of 5.5 at
higher NO3- concentrations (880 uM in f/2 medium). He found that the growth
rate is the same at 100 uM and 880 uM NO3- in the presence of 100 000 to
200 000 cells/mL. Moreover, CO2 limitation does not occur because he finds
the same grpwth rate in cultures bubbled with co2 and cultures agitated by
hand.
He thinks that , as demonstrated by Lomas et
al, at concentrations higher than 100 uM (for instance in f/2 medium with 880 uM NO3-), lots of NO3- is taken up by
diffusive transport, which decrease C:N ratio to 5.5.
Garcia et al (2018) as well as Niemi and Michel (Elementa)found high C:N
ratios (near 7.5) in natural phytoplankton assemblages dominated by
diatoms. In the field, near the SCM, there is around 10 uM NO3- and at the
surface, NO3- concentration decreases more. It thus possible that using a
very rich medium in NO3- induce an artificilly low total C:N ratio because
of luxury consumption of NO3- in vacuoles.

He measured a N:P of 12.5 in the medium with 100 uM NO3- in C. gelidus (u =
0.2 d-1)
Garcia et al measured N:P ratios of 10.10 and a C:N of 5.72 in F. cylindrus
(u = 0.17 d-1)

%}

%{
% If I remove the coupling constants between NADHOR_m and PSI_u, AOX_m = 0,
% but CEF_h = 0.2744

model1 = addReaction(model1, 'PSI_u', 'metaboliteList', {'fdxox_h', 'photon_h', 'pccu1p_u', 'fdxrd_h', 'pccu2p_u', 'mito_coupling_const_h'},...
    'stoichCoeffList', [-2;-2;-2;2;2;-0.15],...  % - 0.15
    'reversible', false);
model1 = addReaction(model1, 'NADHOR_m', 'metaboliteList', {'h_m', 'nadh_m', 'q9_m', 'h_c', 'nad_m', 'q9h2_m', 'mito_coupling_const_h'},...
    'stoichCoeffList', [-5;-1;-1;4;1;1;1],... % 1
    'reversible', false);
model1 = changeRxnBounds(model1, 'AOX_m', 0, 'b');

model1 = addReaction(model1, 'PSI_u', 'metaboliteList', {'fdxox_h', 'photon_h', 'pccu1p_u', 'fdxrd_h', 'pccu2p_u'},...
    'stoichCoeffList', [-2;-2;-2;2;2],...  % - 0.15
    'reversible', false);
model1 = addReaction(model1, 'NADHOR_m', 'metaboliteList', {'h_m', 'nadh_m', 'q9_m', 'h_c', 'nad_m', 'q9h2_m'},...
    'stoichCoeffList', [-5;-1;-1;4;1;1],... % 1
    'reversible', false);
model1 = changeRxnBounds(model1, 'AOX_m', 0, 'l');
model1 = changeRxnBounds(model1, 'AOX_m', 1000, 'u');
%}

%{
Without considering exchange of redox compounds between mitochondria and chloroplast, we found a cyclic electron flow equivalent to a plastoquinone oxidation rate of 0.2744 mmol gDW-1 h-1, which lead to the production of one mol of reduced plastoquinone (pqh2) per mol of substrate or oxidized plastoquinone (pq). This pqh2 is used by cytochrome b6/f complex (CBFC_u), which produces 3.9 protons per mol pqh2. Since 3.97 proton is needed to synthesize one mol ATP by ATPS_h, around 0.98 mol ATP is produced per mol of pq (in CEF_h). Therefore, around 0.27 mmol ATP gDW-1 h-1 is produced through CEF. This represents around 10% of total ATP synthesized in the chloroplast by the ATP synthase. When coupling mitochondrial nadh oxidation and ATP production in the mitochondria to PSI reaction and setting AOX flux equals to 0 (otherwise this reaction is activated and the mitochondrial electron transport chain is short-circuit), the same 0.27 mmol gDW-1 h-1 ATP is produced in the mitochondria by the ATP synthase rather than in the chloroplast, hence simulating efficient communication between both organelles as recently discovered in P. tricornutum.

ATPS_m : 0.239, ATPS_u : 3.24, FNOR_h (2.22851), 


With coupling constant, the model synthesizes more ATP in the mitochondrium using proton gradient generated by the mitochondrial electron transport chain (NADHOR) , which is sustained by AOX_m oxidation of plastoquinone (no fluxes through other electron transport mitochondrial complexes).

No constraints on AOX : ATPS_m : 0.4260 , ATPS_u : 3.384, FNOR_h (2.22851), here excess ATP is produced by ATPS_m and ATPS_u 

If constraints on AOX_m = 0, ATPS_m : 0.52 , ATPS_u : 3.384. Here , there is exactly 0.27-28 mmol ATP produces in the mitochondria to equilibrate ATP:NADPH ratio 

%}



% The N:P ratio is high, but increasing P exchange, makes the problem
% over-constrained in further sensitivity analyses and lead to infeasible
% problems.
%model1 = changeRxnBounds(model1, 'EX_pi_e', -0.004074, 'b');
%model1 = changeRxnBounds(model1, 'EX_no3_e', -0.1383, 'b'); 

% At C:N ratio of 6 and N:P ratio of 28.8, even though net P uptake rate is
% plausible (-0.00444 mmol gDW-1 h-1), gross P uptake rate is very high 0.89 mmol / gDW
% / h through Rxn 2051 (PINAt_e) and gross release is 0.887 mmol / gDW-1
% h-1 through Rxn 1852, PIt_e (h_c + pi_c <=> h_e + pi_e  ). High Cl and h_c efflux
% occurs too through CLt_e (rxn 2049).

% Plasma membrane ATPase (ATPS_c) is inactive. Must add constraints
% This will not lead to pi loss from the cell
% ( h2o_c + atp_c -> h_e + pi_c + adp_c )
% model1 = changeRxnBounds(model1, 'ATPS_c', 0.83, 'l'); % Bd: 0 / +1000
% close pi loss, otherwise, this increase so that Pi loss increase
% model1 = changeRxnBounds(model1, 'PIt_e', 0, 'b'); % Bd: 0 / +1000
% Gives an infeasible problem

% Turning of PEPC_h (phosphoenolpyruvate carboxylase involved in C4
% photosynthesis) since the C4 pathway is incomplete and this consume
% residual CO2.
% model1 = changeRxnBounds(model1,'PEPC_h', 0, 'b'); % Bd: 0 / +1000

% Low mitochondrial respiration at midday (Goldman et al 2015), assume 5%% of all ATP produced : ATPS_h * 5%
% In Fig.4, around 20% (?) of gross photosynthetic O2 evolution is consumed
% in F cylindrus cultures by respiratioon and photorespiratioon.
% Assuming a similar ATP:O2 ratio, which should be around 5 ATP:O2 for
% glucose oxidation and 4.67 ATP:O2 (see rxn ATPS_h) for photosynthesis...
% the relative changes in O2 is comparable to the relative changes in O2
% model1 = changeRxnBounds(model1,'ATPS_m', 1.4, 'b');

% Photorespiration and Mehler reaction small : Goldman et al 2014: assumed
% negligible (optimality assumption)

% Note that The specificity factor (t) of Rubisco, a measure
% of its ability to discriminate CO2 from O2, is considerably higher
% in diatoms than in cyanobacteria and green algae (reviewed in
% [50]) suggesting a lower rate of O2 fixation in diatoms than
% observed in members of the green lineage. This is supported by
% studies showing photorespiratory activity in diatoms at a reduced
% rate than expected from studies with higher plants [51?54]. (see Kroth et
% al 2008 PlosOne).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding constraints on C:N ratio, this does not work. We must change total
% protein concentration instead.

% If I add a constraint on N uptake rate, the model simply choose the
% highest N uptake rate, increase DM_no3_e, but stay constant all other N
% assimilation, if I close 'DM_no3_c', the model excrete nh4 and keep the
% same N:C ratio. If I also close 'EX_nh4_e', the model use other demand
% reaction in the model and grow slowly.
% And if I try constraining directly production of thrtrna_c, cystrna_c, etc., the model
% returns 'INFEASIBLE'.
%N_low =  - Ps / (5.5 * 1.5)  % in mmol / gDW / h % 6.6
%N_high = - Ps / (5.5/1.5)

%model1 = changeRxnBounds(model1,'EX_no3_e', N_high, 'l'); % Bd: -1000 / 0
%model1 = changeRxnBounds(model1,'EX_no3_e', N_low, 'u'); % Bd: -1000 / 0
%model1 = changeRxnBounds(model1,'DM_no3_c', 0, 'l'); % Bd: 0/1000 % BUT NH4 excretion and no increase in N uptake
%model1 = changeRxnBounds(model1,'DM_no3_c', 1000, 'u'); %
%model1 = changeRxnBounds(model1,'EX_nh4_e', 0, 'l');
%model1 = changeRxnBounds(model1,'EX_nh4_e', 1000, 'u');

%model1 = changeRxnBounds(model1,'THRTL_c', FBAsolution.x(264) * 2, 'b');

%FBAsolutionNO3 = optimizeCbModel(model1, 'max', 'one')
%printFluxVector(model1, FBAsolutionNO3.x, false, true)

%{
length(find(FBAsolution.x ~= FBAsolutionNO3.x))
model1.rxns(find(FBAsolution.x ~= FBAsolutionNO3.x))

FBAsolutionNO3.x(findRxnIDs(model1, 'CEF_h')) % 0
FBAsolution.x(findRxnIDs(model1, 'CEF_h')) % 0
FBAsolutionNO3.x(findRxnIDs(model1, 'ATPS_h')) % Same
FBAsolution.x(findRxnIDs(model1, 'ATPS_h')) % Same
FBAsolutionNO3.x(findRxnIDs(model1, 'ATPS_m')) % Same
FBAsolution.x(findRxnIDs(model1, 'ATPS_m')) % Same

find(abs(FBAsolution.x - FBAsolutionNO3.x) > 0.1)

ID = find( abs(FBAsolution.x - FBAsolutionNO3.x) > 1e4 * eps(min(abs(FBAsolution.x),abs(FBAsolutionNO3.x))) )
[num2cell(FBAsolution.x(ID)), num2cell(FBAsolutionNO3.x(ID)), model1.rxns(ID)]
% model1.rxns(ans)
printRxnFormula(model1, model1.rxns(ID));


% becomes active (C4-like mechanisms), Rxn #50  PEPC_h (0 to 0.05074): h2o_h + co2_h + pep_h -> h_h + pi_h + oaa_h 
surfNet(model1, 'PEPC_h', 0, FBAsolutionNO3.x) 
% Rxn #386  PYK_h (0.23285) (decrease from 0.28 to 0.23, by 0.05), h_h + adp_h + pep_h -> atp_h + pyr_h  
surfNet(model1, 'PYK_h', 0, FBAsolutionNO3.x) 
% Rxn #397  PYC_h (0) (inactivated: 0.05 to 0), atp_h + pyr_h + hco3_h ->
% h_h + pi_h + adp_h + oaa_h    TO B13Glucan synthesis...
surfNet(model1, 'PYC_h', 0, FBAsolutionNO3.x) 
%}

%{
% There are a few demand reactions
% DM_no3_c : increases when I increases NO3 uptake. The DM_no3_c also
becomes a higher proportion of NO3 uptake when I increases NO3 uptake.

% A demand reaction for NO3 was added to
account for cellular nitrate that has not yet been assimilated into other biomass components
such as proteins but is included in the dry weight measurements

% DM_phyt_c :  IP52K_c before the demand reaction
% This is the end result of 'Inositol Phosphate metabolism' starting from
% g6p
model.subSystems(findRxnIDs(model, 'IP52K_c'))
SubSystemInositol = ans{1}
Rxns_Inositol = findRxnsFromSubSystem(model, SubSystemInositol)

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_Inositol, 'true', 'struc',...
    {''}, {''})

Rxns_glycolis = findRxnsFromSubSystem(model, SubSystemInositol)

% DM_m2masn_c : MM2_c is the reaction before
% This is the end results of 'N-Glycan biosynthesis' starting with
ipdp/frdp
model.subSystems(findRxnIDs(model, 'MM2_c'))
SubSystemN_Glycan = ans{1}
Rxns_N_Glycan = findRxnsFromSubSystem(model, SubSystemN_Glycan)

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_N_Glycan, 'true', 'struc',...
    {''}, {''})

% DM_thmppp_c, The reaction before is :  TMPPK_c
% So that thmtp_c does not accumulate, a demand reaction is included.

model.subSystems(findRxnIDs(model, 'TMPPK_c'))
SubSystemCofactor_B1 = ans{1}
Rxns_CofactorB1 = findRxnsFromSubSystem(model, SubSystemCofactor_B1)

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_CofactorB1, 'true', 'struc',...
    {''}, {''})

% DM_tre_c, One reaction (TRE6PSP_c) synthesizes 'tre_c'
model.subSystems(findRxnIDs(model, 'TRE6PSP_c'))
SubSystemStarch_Sucr = ans{1}
Rxns_Starch_Sucr = findRxnsFromSubSystem(model, SubSystemStarch_Sucr)

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_Starch_Sucr, 'true', 'struc',...
    {''}, {''})
% This small subsystem synthesises UDP-glucuronate, which is used for the carbohydrate biomass synthesis
% and trehalose, which leaves the system through a demand reaction.

% DM_5drib_c : 5DOAN_c produces 5dribc (5'-deoxyribose), which leaves the
system.
model.subSystems(findRxnIDs(model, '5DOAN_c'))
SubSystemBiotin = ans{1}
Rxns_Biotin = findRxnsFromSubSystem(model, SubSystemBiotin)

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_Biotin, 'true', 'struc',...
    {''}, {''})
% dad_5 (deoxyadenosine) is produced when biotin is synthesized in the mitochodnria and when
% thiamin is synthesized in the chloroplast. Then it is converted in 5drib,
which leaves the system

% DM_2mop_m : 2mop_m can produces ppcoa, which is then convert in organic
acid, succoa for the TCA cycle. However, I think that 2mop_m can accumulate if degration
is present , but the TCA cycle cannot use all incoming organic acid.

% DM_fald_m : fald_m (formaldehyde) is produced when sarcosine is produced,
but no loss reaction is present. fald_m must leave the system.

% sink_Asn-X-Ser_Thr_c : The reaction 'DOLASNT_c,' uses the species Asn-X-Ser_Thr_c
% This produces Asn-X-Ser_Thr_c from nothing if needed... This reaction
should be equal to 0....
model.subSystems(findRxnIDs(model, 'DOLASNT_c'))
SubSystemN_Glycan = ans{1}
Rxns_N_Glycan = findRxnsFromSubSystem(model, SubSystemN_Glycan)

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_N_Glycan, 'true', 'struc',...
    {''}, {''})

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_Folate, 'true', 'struc',...
    {''}, {''})

% DM_mhpglu_c : 5mthglu_c is produced by MTHF3ES_c
model.subSystems(findRxnIDs(model, 'MTHF3ES_c'))
SubSystem_Folate = ans{1}
Rxns_Folate = findRxnsFromSubSystem(model, SubSystem_Folate)

[involvedMets, deadEnds] = draw_by_rxn(model, Rxns_Folate, 'true', 'struc',...
    {''}, {''})
% mthf produces 5mthf, which produces 5mthglu_c, and leaves the system

% DM_dmsp_c
% Allow growth on cysteine, which could occur by a putative conversion of
% DMSP into amino acids.

% DM_indole_c 
% Indole accumulation prohibited growth on tryptophan as nitrogen source. To account
% for the unknown indole degradation, a demand reaction was added

% DM_for_c : Allow histidine catabolism
% Since we could not identify genes
% that are involved in histidine catabolism in P. tricornutum, we added histidine catabolism as
% one lumped, low confidence reaction degrading histidine and water into ammonium, formamide
% and glutamate. Formamide is split into formate and ammonium with formate accumulating
% during histidine catabolism in silico; a demand reaction was added to allow the accumulated
% formate to leave the system.

%}


%%%%%%%%%%%%%%%%%%%%

% C uptake rate = 0.5 g C / g DW / 0.3 d = 1.66 g C / g DW / d = 0.0692 g
%  C / g DW / h or 0.0058 mol C g DW / h

% 0.5 g C / g DW * 0.3 d-1 = 0.15 g C / g DW / d = 0.0062 g C / gDW
% / h or 0.5167 mmol C / g DW / h

% A Redfield C:N:P ratio is assumed (106:16:1)
% The HFSP measured C: N ratio are close to 6.6 (106:16)
% But, can be close to 5.5 - 6 in polar diatoms...
% CO2_fixed = -((0.5 / 12 * 1000) * 0.3 / 24) % in mmol / gDW / h






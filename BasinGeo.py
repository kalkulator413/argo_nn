import shapely
import matplotlib.pyplot as plt

# taken from Paul's code here: https://github.com/Chamberpain/OptimalArray/blob/main/Utilities/CorGeo.py

basins = {
	'Indian':[(147,-43.52940220185144),(143.3479485275823,-37.33601703562596),(129.4747260606591,-31.65737138513975),
        (116.8213348808513,-33.63996335620772),(114.1336886317196,-24.11210578182048),(126.8895584292698,-20.10475922355295),
        (116.1158834870691,-5.730689136449962),(106.3794488756852,-6.953056836849058),(96.28600509943337,4.831929831532299),
        (97.90360579962756,16.87704560343002),(94.27317289741688,16.37361837248965),(94.33359097892006,19.37403360157845),
        (90.2091678810215,21.29516387034062),(85.33933054916371,20.28690337126939),(79.90939997726753,15.28725891185891),
        (79.73043452807434,12.76560688789527),(77.61698366176469,8.677663553177529),(73.07173102493962,17.82552817022485),
        (72.45900019949667,21.53417968653654),(69.83928025251217,21.52699330766178),(67.39703210720685,24.5291555329574),
        (66.46821530686643,25.78691637280379),(56.70750442437404,25.66075493384507),(58.72984841152353,22.03163388871144),
        (57.06581416303475,19.15932876978664),(52.34927042817224,16.21465648495196),(43.4652228713137,12.71268431010668),
        (41.93132280452403,10.95937698688315),(50.01226907434041,9.405096543084687),(46.12570918863462,3.933998362277269),
        (41.20631927073987,-0.9345519525315622),(37.60478119540351,-6.210139259824268),(39.67592988828236,-15.07135869457443),
        (33.69621689976939,-19.65087962850418),(34.1381515444963,-23.93162439597786),(30.02506496520469,-29.23313789495232),
        (26.55519092843285,-33.21902230207225),(20.0,-34.34985699240418)],
	'NAtlantic':[(-15.99298965537844,20.0),( -16.51768345340904,21.5135791895485),( -12.85930824574853,27.3984923294593),
        ( -9.814281489721248,27.97722929568377),( -8.76926137454598,31.42705281354388),( -4.844078566635137,33.22829780277818),
        ( -7.163320253103191,43.50599870860887),( -1.036215450782537,42.18863705314655),( 0.9480128425790557,43.98587477724783),
        ( -2.712196087815386,47.86200410331944),( -9.490699205550818,52.24203463882663),( -1.430837126611296,60.492833684356),
        ( -18.48518831195036,64.94737597490287),( -26.693204331776,69.74009059144014),( -46.43161024784422,61.27338745482746),
        ( -50.20187747100773,72.80465043844811),( -67.83861820446973,78.25341280867181),( -84.00411784512445,74.95668902978969),
        ( -63.88158735865454,66.14474917409777),( -64.83285124543914,60.15569777481804),( -52.62163270584416,48.60704311039531),
        ( -61.17606612066979,45.85551735135405),( -74.15449943541451,41.02571986947127),( -75.97004374470255,38.21535215461827),
        ( -75.74775841682803,35.53781825873987),( -80.78884359192395,32.69703604266579),( -81.86666651272962,30.71636862921097),
        ( -80.5420695665382,25.70069466689458),( -80.7262399142712,22.75689399043104),( -75.62705224575633,20.64254575052694),
        ( -72.78563669764426,20.0),( -15.99298965537844,20.0)],
	'TropicalAtlantic':[(-72.50562296593817,20),(-70.59819414971815,18.42528248241784),(-64.77041104794849,18.32119324306983),
        (-61.84634301084144,17.50229847613128),(-59.59116865706082,14.6807045454499),(-59.50143599673747,11.99113053869993),
        (-61.83274263008427,9.136119519346721),(-57.27276214022088,5.500062512667464),(-53.46534316823615,5.052030577588845),
        (-51.58589087805466,3.537188235266167),(-50.43166199124011,0.9209177789148023),(-48.9658815572249,-0.9123304819710399),
        (-42.79308489712654,-3.664005213856391),(-39.96737056590514,-3.510326760772562),(-35.99174323989286,-6.158377106084419),
        (-39.30672759560186,-11.90655999832168),(-41.69766615837543,-20.0),(13.09374881888988,-20.0),
        (12.11794803335862,-17.08155116752772),(14.31502302784277,-11.29020044115646),(13.64063555778271,-7.507139436739451),
        (9.671662939732357,-1.228970263903303),(10.8276175212515,3.634534464349576),(5.210162815943946,6.789548635201794),
        (2.266706807472187,6.34039712947654),(-4.43823577108877,5.292657896690207),(-7.670451218240339,3.969893556481763),
        (-12.87968607914255,7.743606278363659),(-16.4467533222392,12.22832083758344),(-16.99784648828647,14.6786385909283),
        (-15.71191486006144,18.04994037837994),(-15.97131585254839,20.0),(-72.50562296593817,20)],
	'SouthAtlantic':[(20.0,-34.35046375144582),(18.67347013265745,-33.96063336022314),(18.02848113780506,-32.94443890987839),
		(18.49489620471581,-32.56463790467472),(18.33072742274843,-31.72700300431016),(16.89790192337395,-28.91744290987718),
		(15.79225162404466,-27.98101912715219),(15.32499178636868,-26.87953225177936),(15.06980111871543,-26.27891323981478),
		(14.95293298933232,-25.02540654629124),(14.68418587154647,-24.53932913343425),(14.5462570131453,-23.03528728932478),
		(14.47892381265254,-22.33408187910103),(14.04831979105372,-21.77986345796982),(13.43194958268201,-20.82128066772663),
		(13.28859012755361,-20.33722934812156),(13.09363193440577,-20.0000),(-40.5976304911968,-20.0),
		(-42.06014849903111,-22.18865116561715),(-49.08751869488754,-24.61704752875797),(-49.09437772573018,-27.71168486619063),
		(-53.95184717028173,-33.76925336898484),(-59.32745675268787,-33.78030719803291),(-57.57852470583949,-37.29365804156026),
		(-62.55109245211934,-38.43858253207481),(-62.76054572116809,-40.57152851528756),(-65.64245834496685,-40.1785373343289),
		(-65.59474035785468,-44.22600517770126),(-68.23921716618344,-46.07976368124183),(-66.33445934267897,-47.61800182751628),
		(-69.44497824762939,-50.47718723831122),(-69.17545169568312,-53.05305794039962),(-66.00425784471349,-54.75603363189224),
		(-67.0,-56.03247111630986)],
    'NPacific':[(-130.4035261964233, 55.0), (-127.0812674058722, 49.62755804643379), (-123.5031152865783, 47.69525998053675),
        (-123.9821609070654, 42.58507262179968), (-123.4646493631059, 38.58314108287027), (-121.8138182188833, 35.72586240471471),
        (-117.3894585799435, 32.49679570904591), (-113.2985172073168, 31.65136326179817), (-105.4266560754428, 23.05901404803846),
        (-104.6431889409656, 20.0), (-135.0, 20.0), (-180.1, 20.0), (-180.1, 51.8), (-171.7988217450257, 52.44580890664249),
        (-155.5893394435851, 59.05068402489769), (-146.447784048455, 61.83947608617396), (-130.4464882993783, 55.76859441423366),
        (-130.5722813264113, 55.0), (-130.4035261964233, 55.0)],
    'TropicalPacific':[(-180.1, -20), (-68.42930684879389, -20), (-70.17015104761151, -16.55216215970295),
        (-74.96754093664958, -14.03860742174251), (-79.73669189883009, -5.584892759052865), (-78.52106685157878, 0.408381155896844), 
        (-76.76810851638021, 7.663643325758725), (-79.35212462962738, 9.1729278834978), (-82.45222772097509, 8.188702540343506), 
        (-84.42711652089258, 10.32739929197514), (-88.88605864053, 14.11861449658679), (-93.83727342627287, 16.40515817053415), 
        (-104.28851657329369, 20), (-180.1, 20), (-180.1, -20)],
    'SouthPacific':[(-67, -56.05831143635918), (-70.01906226899843, -54.94099908902054), (-74.4700954236148, -50.6355821528286),
        (-72.95976403442722, -42.09748969048012), (-70.60667601995038, -32.89435779688432), (-68.99964469056277, -20.0), (-180.1, -20)],
    'GOM': [(-93.83443543373453,30.31192409716139),(-97.33301191109445,28.53461995311981),(-98.22951964586647,26.56102839784981),
        (-98.23176621309408,21.14640430279481),(-94.75173730705114,17.69018757346009),(-89.11465568491809,17.20328330857444),
        (-89.06712280579349,15.47270523764778),(-84.25124860521271,15.32356939252626),(-84.11167583129146,10.60002449442026),
        (-81.68974915579967,8.513630710186039),(-79.21737886089964,9.102618171647688),(-76.35489819907755,7.748300845902937),
        (-72.52740327179197,10.70167332395043),(-62.02934864430507,9.44397718074659),(-59.67683898413669,12.02173132700068),
        (-59.82535447869236,14.70093873644528),(-61.79537355504573,17.49727742725868),(-64.89185977095053,18.25534196301828),
        (-70.61882162015974,18.45320848406356),(-72.82039923208698,19.93626561168838),(-80.81793590533317,22.75362929988757),
        (-80.70017685693531,25.6944945895363),(-82.68021005951927,30.09720890646617),(-93.83443543373453,30.31192409716139)]
}

# class InverseSO(InverseGeo):
# 	facecolor = 'plum'
# 	facename = 'Southern Ocean'
# 	plot_class = SOSECartopy
# 	region = 'southern_ocean'
# 	lat_sep=1
# 	lon_sep=1
# 	l=3

# 	def __init__(self,*args,**kwargs):
# 		file = os.path.join(DATA_DIR,'saf.asc')
# 		token = open(file,'r')
# 		coord = [i.strip().split() for i in token.readlines()]
# 		xcoord,ycoord = zip(*[(float(i[0]),float(i[1])) for i in coord if abs(float(i[0]))<179])	
# 		xarray = np.array(xcoord)
# 		yarray = np.array(ycoord)
# 		coord_list = list(zip(xcoord, ycoord))
# 		coord_list += [(180,ycoord[-1]),(180,-90),(-180,-90),(-180,ycoord[0]),coord_list[0]]
# 		self.coord_list = coord_list
# 		super().__init__(*args,**kwargs)

if __name__ == '__main__':
    for key in basins:
        x,y = shapely.Polygon(basins[key]).exterior.xy
        plt.plot(x,y)
    # plt.legend(basins.keys())
    plt.show()

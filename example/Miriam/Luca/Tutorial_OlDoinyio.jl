using DelimitedFiles, GeophysicalModelGenerator, JLD2
p               = ProjectionPoint(Lon=35.76, Lat=-3.01);

# 0. load topography and geological map
Topo_Degrees    = load("Topo_Oldoinyo.jld2","Topo_midres");
Lon             = Topo_Degrees.lon.val;
Lat             = Topo_Degrees.lat.val;
Depth           = Topo_Degrees.depth.val;
T               = Topo_Degrees.fields.Topography;
Topo_Degrees    = GeophysicalModelGenerator.GeoData(Lon, Lat, T*1000, (Topography = T,));
Topo_Cart       = Convert2CartData(Topo_Degrees,p);
Write_Paraview(Topo_Degrees,"Ol_Topo_Degrees")
Write_Paraview(Topo_Cart,"Ol_Topo_Cart")

Corner_LowerLeft   = ( 35.76, -3.01, 0)
Corner_UpperRight  = ( 36.29, -2.50, 0)
mappa           = Screenshot_To_GeoData("cut_geology.png", Corner_LowerLeft, Corner_UpperRight)
mappa_Cart      = Convert2CartData(mappa,p);
Write_Paraview(mappa, "geological map")
Write_Paraview(mappa_Cart, "geologicalCart")

# 1. load earthquake data and save in degrees
data                = readdlm("locations.txt", ',', Float64, skipstart=0, header=false);
lat                 = data[:,1];
lon                 = data[:,2];
depth               = -data[:,3]/1000;
EQ_Data_Degrees     = GeophysicalModelGenerator.GeoData(lon, lat, depth, (Depth = -depth/1000 * km,));
Write_Paraview(EQ_Data_Degrees, "Ol_Earthquakes_Degrees", PointsData=true);
depth               = -data[:,3];
EQ_Data_Degrees     = GeophysicalModelGenerator.GeoData(lon, lat, depth, (Depth = -depth/1000 * km,));
EQ_Data_Cart        = Convert2CartData(EQ_Data_Degrees,p);
Write_Paraview(EQ_Data_Cart, "Ol_Earthquakes_Cart", PointsData=true);

# 2. load the attenuation models and save
data_pd             =   readdlm("peakdelay_6_Degrees_Hz.txt", ',', Float64, skipstart=0, header=false);
data                =   readdlm("Qc6Hz.txt", ',', Float64, skipstart=0, header=false);
lon                 =   data[:,1];
lat                 =   data[:,2];
depth               =   data[:,3]/1000;
peakDelay           =   data_pd[:,4];
Qc                  =   data[:,4];
resolution          =   (length(unique(depth)),  length(unique(lat)), length(unique(lon)))
dim_perm            =   [3 2 1]
Lon                 =   permutedims(reshape(lon, resolution), dim_perm);
Lat                 =   permutedims(reshape(lat, resolution), dim_perm);
Depth               =   permutedims(reshape(depth, resolution), dim_perm);
p3d                 =   permutedims(reshape(peakDelay, resolution), dim_perm);
Q3d                 =   permutedims(reshape(Qc, resolution), dim_perm);
Data_set_Degrees    =   GeophysicalModelGenerator.GeoData(Lon, Lat, Depth, (Qc = Q3d, PeakDelay = p3d))
Write_Paraview(Data_set_Degrees, "Ol_Qc_Degrees")
Depth               =   Depth*1000;
Data_set_Degrees    =   GeophysicalModelGenerator.GeoData(Lon, Lat, Depth, (Qc = Q3d, PeakDelay = p3d))
Data_set_Cart       =   Convert2CartData(Data_set_Degrees,p);
Write_Paraview(Data_set_Cart, "Ol_Qc_Cart")

# 3. load the velcity models and save
data_v              =   readdlm("Oldoinyo_modVp.txt", ',', Float64, skipstart=0, header=false);
lon                 =   data_v[:,1];
lat                 =   data_v[:,2];
depth               =   data_v[:,3];
v                   =   data_v[:,4];
resolution          =   (length(unique(depth)),  length(unique(lat)), length(unique(lon)))
dim_perm            =   [3 2 1]
Lon                 =   permutedims(reshape(lon, resolution), dim_perm);
Lat                 =   permutedims(reshape(lat, resolution), dim_perm);
Depth               =   permutedims(reshape(depth, resolution), dim_perm);
V                   =   permutedims(reshape(v, resolution), dim_perm);
Data_set_Degrees    =   GeophysicalModelGenerator.GeoData(Lon, Lat, Depth, (Velocity = V,))
Write_Paraview(Data_set_Degrees, "Ol_Vp_Degrees")
Depth               =   Depth*1000;
Data_set_Degrees    =   GeophysicalModelGenerator.GeoData(Lon, Lat, Depth, (Velocity = V,))
Data_set_Cart       =   Convert2CartData(Data_set_Degrees,p);
Write_Paraview(Data_set_Cart, "Ol_Vp_Cart")


# # Seismicity movie - for now, COMMENTED
# p2                  = @__FILE__;
# p2last              = findlast("/",p2);
# p3                  = chop(p2,head=0,tail = length(p2)-p2last[1]+1);
# output_path         = string(p3,"/");
# movie_p             = Movie_Paraview(name=string(p3,"/TemporalSeismicity"), Initialize=true);
# if isdir(string(p3,"/TemporalSeismicity"))==0
#     mkdir(string(p3,"/TemporalSeismicity"));
# end

# nt                  = 100;
# ldf                 = length(depth)
# for itime = 1:nt:ldf-nt
#     name            = string(p3,"/TemporalSeismicity/", string(itime));
#     la              = lat[itime:itime+nt-1];
#     lo              = lon[itime:itime+nt-1];
#     de              = depth[itime:itime+nt-1];
#     label_time      = itime;
#     Data_set        = GeophysicalModelGenerator.GeoData(lo, la, de, (Depth = de * km,));
#     movie_p         = Write_Paraview(Data_set, name, pvd=movie_p, time=label_time, PointsData=true);
    
# end
# Movie_Paraview(pvd=movie_p, Finalize=true)

import gnss_lib_py as glp

sp3_path = "data/COD0MGXFIN_20211180000_01D_05M_ORB.SP3"
sp3 = glp.Sp3(sp3_path)
print(sp3)

sp3_first_ten_gps = sp3.where("gnss_id","gps")
fig = glp.plot_metric_by_constellation(sp3_first_ten_gps,"gps_millis","x_sv_m")
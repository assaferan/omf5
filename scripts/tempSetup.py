from lmfdb import db
os.chdir("../omf5_data/hecke_evs_3_0/data")
gustavo = eval(open("../../newforms30_dims.txt").read())
gustavo = [[]] + gustavo
os.chdir("../../../omf5")
from scripts.parse_omf import load_ev_data
os.chdir("../omf5_data/hecke_evs_3_0/data")
ev_data_100, G_type_100 = load_ev_data(3,0,"./")
all_dims_100 = [sorted([len(f['field_poly'])-1 for f in ev_data_100[N] if f['aut_rep_type'] == 'G']) for N in range(1000)]
from scripts.parse_omf import update_ev_data_class
update_ev_data_class(ev_data_100, 1000, db)


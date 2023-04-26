def is_yoshida(q, alpha_q_2, alpha_q_4, ap_array, N):

    res2 = db.mf_newforms.search({"level" : {"$in" : divisors(N)}, "weight" : 2, "char_order" : 1}, "hecke_orbit_code")
    res4 = db.mf_newforms.search({"level" : {"$in" : divisors(N)}, "weight" : 4, "char_order" : 1}, "hecke_orbit_code")

    orbits2 = [hoc for hoc in res2]
    orbits4 = [hoc for hoc in res4]
    
    res2_q = db.mf_hecke_traces.search({"hecke_orbit_code" : {"$in" : orbits2}, "n" : q, "trace_an" : alpha_q_2}, "hecke_orbit_code")
    res4_q = db.mf_hecke_traces.search({"hecke_orbit_code" : {"$in" : orbits4}, "n" : q, "trace_an" : alpha_q_4}, "hecke_orbit_code")
    
    orbits4_q = [hoc for hoc in res4_q]
    orbits2_q = [hoc for hoc in res2_q]
    
    ap_4 = [{a["n"] : a["trace_an"] for a in db.mf_hecke_traces.search({"hecke_orbit_code" : hoc4})} for hoc4 in orbits4_q]
    ap_2 = [{a["n"] : a["trace_an"] for a in db.mf_hecke_traces.search({"hecke_orbit_code" : hoc2})} for hoc2 in orbits2_q]
    lamda_ps = [[p*aap_2[p] + aap_4[p] for p in prime_range(100) if N % p != 0] for aap_2 in ap_2 for aap_4 in ap_4]
    return any([lamda_p == ap_array for lamda_p in lamda_ps])

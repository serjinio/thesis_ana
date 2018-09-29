CXX = g++
CPPFLAGS += `root-config --cflags --glibs` -I/Users/serj/ana/include

LDLIBS   += `root-config --glibs` -lMathMore -lProof -lTree -lTreePlayer

ALL: espri_nais_ext espri_nais_double_hit_evts espri_abang_compare rdc_eff2

asymmetry_he4_wo_bg: asymmetry_he4_wo_bg.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $^ -o $@ $(LDLIBS)

espri_nais_ext: espri_nais_ext.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $^ -o $@ $(LDLIBS)

espri_nais_double_hit_evts: espri_nais_double_hit_evts.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $^ -o $@ $(LDLIBS)

rdc_eff2: rdc_eff2.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $^ -o $@ $(LDLIBS)

espri_abang_compare: espri_abang_compare.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $^ -o $@ $(LDLIBS)

clean:
	rm -f espri_e_calib_para asymmetry_he4_wo_bg espri_nais_ext espri_nais_double_hit_evts espri_abang_compare rdc_eff2 *.o

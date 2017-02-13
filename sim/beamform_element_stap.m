function [w v] = beamform_element_stap(D,Vs,Npulses,snapshots,guard_bandwidth,points_kept)
[w v] = element_stap_2(D,Vs,Npulses,snapshots,guard_bandwidth,points_kept);
Quick Help:

To run VLM_wp.m simply download all the files, ensure they are in a single folder.  
Hit run in the editor tab.

To use as a function:

Ensure that the files: add_angle_of_attack.m, add_dihedral.m, add_sweep.m, bf2wf.m, chord_x.m, and horse.m are in the same folder as VLM_f.m.
That way if you want to jointly use VLM_f.m with other code VLM_f.m has the support functions with it.  Simply use:

[CL, CDi, e] = VLM_f(AR, lambda, N, sweep, alfa, beta, Uinf, errors_on)

as the function call.
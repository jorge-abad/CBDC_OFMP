function ds = dynamic_set_auxiliary_series(ds, params)
%
% Computes auxiliary variables of the dynamic model
%
ds.AUX_DIFF_FWRD_11=ds.pi-ds.pi(-1);
ds.AUX_DIFF_FWRD_4=ds.C-ds.C(-1);
ds.AUX_DIFF_FWRD_2=ds.Q-ds.Q(-1);
ds.AUX_DIFF_FWRD_3=ds.I-ds.I(-1);
ds.AUX_DIFF_FWRD_13=ds.Xi1-ds.Xi1(-1);
ds.AUX_DIFF_FWRD_14=ds.Xi2-ds.Xi2(-1);
ds.AUX_DIFF_FWRD_16=ds.Ra-ds.Ra(-1);
ds.AUX_DIFF_FWRD_43=ds.V-ds.V(-1);
end

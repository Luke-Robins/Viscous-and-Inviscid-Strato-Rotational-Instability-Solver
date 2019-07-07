function [ E_, fig_vertical, fig_horizontal ] = SRI_visualise( eta_, mu_, m_, k_, Re_, N_, T_ )
%Dr Luke Robins 2019 luke.robins@cantab.net
%

[ E_, u_k, v_k, w_k, ~, ~ ] = SRI_solver( eta_, mu_, m_, k_, Re_, N_, T_ );

[ fig_vertical, fig_horizontal ] = Visualise( u_k, v_k, w_k, eta_, mu_, m_, k_, N_ );

end
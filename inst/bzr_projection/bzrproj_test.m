%% BZRPROJ_TEST
%   From "Bézier projection: A unified approach for local
%   projection and quadrature-free refinement and coarsening of NURBS and
%   T-splines with particular application to isogeometric design and
%   analysis", D.C. Thomas et al., Comput. Methods Appl. Mech. Engrg. 284
%   (2015)55-105
%
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

clc;
clear;
close all

%% Refine curve------------------------------------------------------------
% curve parameters
p = 2;
controlPoints = [[-0.1 -0.3 -0.15 0.7 0.9 1];[0.2 0.4 1.0 0.9 0.75 0.4]];
source_knots = [0 0 0 1/4 1/2 3/4 1 1 1];

% plot
crv = nrbmak(controlPoints, source_knots);
figure(1); 
nrbplot(crv, 200);
hold on;
grid on;

%% Bézier projection element level-----------------------------------------
% construct the Bézier extraction operator for the source and the target
% mesh

target_knots = [0 0 0 1/4 1/2 3/4 7/8 1 1 1];
refined_knots = [0 0 0 1/8 1/4 3/8 1/2 5/8 3/4 7/8 1 1 1];
C = bzrextr(source_knots, p);
C_coarse = bzrextr(target_knots, p);
C_refined = bzrextr(refined_knots, p);
complete_contolpoints = cat(2, [0 0 0 0 0 0 0 0 0 0]', [1 1 1 1 1 1 1 1 1 1]');

Bref_1 = bzrproj_el_href( p, C(:,:,1), C_refined, [1], 2 );
Bref_2 = bzrproj_el_href( p, C(:,:,1), C_refined, [2], 2 );
Bref_3 = bzrproj_el_href( p, C(:,:,2), C_refined, [3], 2 );
Bref_4 = bzrproj_el_href( p, C(:,:,2), C_refined, [4], 2 );
Bref_5 = bzrproj_el_href( p, C(:,:,3), C_refined, [5], 2 );
Bref_6 = bzrproj_el_href( p, C(:,:,3), C_refined, [6], 2 );
Bref_7 = bzrproj_el_href( p, C(:,:,4), C_refined, [7], 2 );
Bref_8 = bzrproj_el_href( p, C(:,:,4), C_refined, [8], 2 );

ref_control_point_e1 = Bref_1(:,:,1) * crv.coefs(1:2, 1:3)';
ref_control_point_e2 = Bref_2(:,:,1) * crv.coefs(1:2, 1:3)';
ref_control_point_e3 = Bref_3(:,:,1) * crv.coefs(1:2, 2:4)';
ref_control_point_e4 = Bref_4(:,:,1) * crv.coefs(1:2, 2:4)';
ref_control_point_e5 = Bref_5(:,:,1) * crv.coefs(1:2, 3:5)';
ref_control_point_e6 = Bref_6(:,:,1) * crv.coefs(1:2, 3:5)';
ref_control_point_e7 = Bref_7(:,:,1) * crv.coefs(1:2, 4:6)';
ref_control_point_e8 = Bref_8(:,:,2) * crv.coefs(1:2, 4:6)';

ref_control_points(1,:) = ref_control_point_e1(1,:);
ref_control_points(2,:) = ref_control_point_e1(2,:);
ref_control_points(3,:) = ref_control_point_e1(3,:);
ref_control_points(4,:) = ref_control_point_e3(2,:);
ref_control_points(5,:) = ref_control_point_e3(3,:);
ref_control_points(6,:) = ref_control_point_e5(2,:);
ref_control_points(7,:) = ref_control_point_e5(3,:);
ref_control_points(8,:) = ref_control_point_e7(2,:);
ref_control_points(9,:) = ref_control_point_e7(3,:);
ref_control_points(10,:) = ref_control_point_e8(3,:);

refcrv = nrbmak((cat(2, ref_control_points, complete_contolpoints))', refined_knots);
nrbplot(refcrv, 200);

% project control values of the source mesh onto the target mesh
B_1 =  bzrproj_el_hcoarse( p, C_refined, inv(C_coarse(:,:,1)), 1/8, 1/4, [1 2] );
B_2 =  bzrproj_el_hcoarse( p, C_refined, inv(C_coarse(:,:,2)), 1/8, 1/4, [3 4] );
B_3 =  bzrproj_el_hcoarse( p, C_refined, inv(C_coarse(:,:,3)), 1/8, 1/4, [5 6] );
B_4 =  bzrproj_el_hcoarse( p, C_refined, inv(C_coarse(:,:,4)), 1/8, 1/8, 7 );
B_5 =  bzrproj_el_hcoarse( p, C_refined, inv(C_coarse(:,:,5)), 1/8, 1/8, 8 );

complete_contolpoints = cat(2, [0 0 0 0 0 0 0]', [1 1 1 1 1 1 1]');

% new Bézier elements control points
new_control_point_e1 = B_1(:,:,1) * ref_control_points(1:3, 1:2) + B_1(:,:,2) * ref_control_points(2:4, 1:2);
new_control_point_e2 = B_2(:,:,1) * ref_control_points(3:5, 1:2) + B_2(:,:,2) * ref_control_points(4:6, 1:2);
new_control_point_e3 = B_3(:,:,1) * ref_control_points(5:7, 1:2) + B_3(:,:,2) * ref_control_points(6:8, 1:2);
new_control_point_e4 = B_4(:,:,1) * ref_control_points(7:9, 1:2);
new_control_point_e5 = B_5(:,:,1) * ref_control_points(8:10, 1:2);


% smoothing using local least square weights
new_control_points(1,:) = new_control_point_e1(1,:);
new_control_points(2,:) = 1/2 * new_control_point_e1(2,:) + 1/2 * new_control_point_e2(1,:);
new_control_points(3,:) = 1/3 * new_control_point_e1(3,:) + 1/3 * new_control_point_e2(2,:) + 1/3 * new_control_point_e3(1,:);
new_control_points(4,:) = 1/3 * new_control_point_e2(3,:) + 1/3 * new_control_point_e3(2,:) + 1/3 * new_control_point_e4(1,:);
new_control_points(5,:) = 1/3 * new_control_point_e3(3,:) + 1/3 * new_control_point_e4(2,:) + 1/3 * new_control_point_e5(1,:);
new_control_points(6,:) = 1/2 * new_control_point_e4(3,:) + 1/2 * new_control_point_e5(2,:);
new_control_points(7,:) = new_control_point_e5(3,:);


% make B-spline
pcrv = nrbmak((cat(2, new_control_points, complete_contolpoints))', target_knots);

nrbplot(pcrv, 200);
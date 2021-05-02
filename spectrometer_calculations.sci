// Copyright (C) 2021 - Buğra Coşkun
// bgra.coskn@gmail.com  ||  https://bgracoskn.github.io
//
// This work is licensed under the Creative Commons 
// Attribution-ShareAlike 4.0 International License. 
// To view a copy of this license, visit 
// http://creativecommons.org/licenses/by-sa/4.0/ 
// or send a letter to Creative Commons, PO Box 1866, 
// Mountain View, CA 94042, USA.
//
// Date of creation: 2021
//
// Related links:
// Instructables, github, thingiverse, bgracoskn, linkedin
clc; clear all; funcprot(0); format(25)
//atomsLoad('DD_QD')
wl_redlimit = 0.7; //700nm violet
wl_violetlimit = 0.380; //380nm red

function deg = rad2deg(rads)
	deg = (360./(2.*%pi)).*rads;
endfunction

function rad = deg2rad(degs)
	rad = (%pi/180)*degs
endfunction

function na = ref_ind_air(wla)
	//https://refractiveindex.info/?shelf=other&book=air&page=Borzsonyi
	part1 = (14926.44*(10^(-8))*(wla^2))/(wla^2-19.36*(10^(-6)));
	part2 = (41807.57*(10^(-8))*(wla^2))/(wla^2-7.434*(10^(-3)));
	na = sqrt(part1 + part2 + 1)
endfunction

function np = ref_ind_prism(wlp)
	//https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Szczurowski
	part1 = (0.99654*wlp^2)/(wlp^2-0.00787);
	part2 = (0.18964*wlp^2)/(wlp^2-0.02191);
	part3 = (0.00411*wlp^2)/(wlp^2-3.85727);
	np = sqrt(part1+part2+part3+1)
endfunction

function _paths_ = p_path(phi, al, na, np)
	// Isosceles triangle prism path calculation function.
	// deg(OAB) = deg(OBA) = beta ==> deg(AOB) = 180 - 2beta = al
	// phi = OA orthogonal surface angle to i or x axis
	// d0 = incoming light ray angle, if parallel to x: d0=phi
	// al: assuming isosceles triangle, al is the angle between equal sides
	// _d2_ is the angle of leaving light ray with the i or x axis
	d0 = %pi/2 - phi
	d1 = asin((na*sin(d0))/np)
	d1_ = al - d1
	d2 = asin((np*sin(d1_))/na) // angle between surface normal and ray
	_d2_ = %pi/2 - phi - al - d2 // angle of ray to horizontal axis
	disp(_d2_)
	_paths_ = [d0; d1; d1_; d2; _d2_]
endfunction

function _output_ = cal_path(phi,al,wlr,wlb, h1)
	//h1 is the distance of ray in the vertical axis to the top point of prism
	//h2 is the distance of leaving ray in the vertical axis to the top point of prism
	// $\theta_0 =\dfrac{\pi}{2} - \phi \\ n_a\sin(\theta_0)=n_p\sin(\theta_1) \\ \theta_1^'=\alpha-\theta_1 \\ n_p\sin(\theta_1^')=n_a\sin(\theta_2) \\ \theta_2^'=\pi/2-\phi-\alpha+\theta_2$
	// $ h_2=\dfrac{h_1 \cos(\theta_1^') \sin(\alpha + \phi)}{\cos(\theta_1)\sin(\phi)}}$
	na_red = ref_ind_air(wlr); 
	disp('na_red', na_red)

	na_blue = ref_ind_air(wlb); 
	disp('na_blue', na_blue)

	np_red = ref_ind_prism(wlr); 
	disp('np_red', np_red)

	np_blue = ref_ind_prism(wlb); 
	disp('np_blue', np_blue)

	red_path = p_path(phi, al, na_red, np_red);
	red_h2 = h1*cos(red_path(3))*sin(al+phi) / ..
		(sin(phi)*cos(red_path(2)))
	
	blue_path = p_path(phi, al, na_blue, np_blue);
	blue_h2 = h1*cos(blue_path(3))*sin(al+phi) / ..
		(sin(phi)*cos(blue_path(2)))
	
	_output_ = [red_path(5);blue_path(5);red_h2;blue_h2]
	disp('output:',_output_)
	
	disp('red_path:', rad2deg(red_path))
	disp('blue_path:', rad2deg(blue_path))
	
endfunction

output_vals = cal_path(deg2rad(67.5), %pi/4, wl_redlimit, wl_violetlimit, 10) //h1=10mm 
output_vals_deg = [rad2deg(output_vals(1));rad2deg(output_vals(2));output_vals(3:4)]
disp(output_vals_deg)

// print a testing rig before going to the concave mirror calculations









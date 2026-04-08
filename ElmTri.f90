program ElmTri 

! VERSION WITH TRIANGULAR END SLICES AND NORMAL FORCE AT DX/3
! Fortran F90 program for the calculation of the forces on a single
! fully specified potential slip surface by the enhanced limit method
! comprised of triangular slices in accordance with the method
! proposed by Alexandre and Thomasi (2025)

! Compilation in terminal linux:
! gfortran ElmTri.f90 -o ElmTri

! Compilation in terminal linux with debugging steps:
! gfortran ElmTri.f90 -O0 -fcheck=all -Wall -Wextra -fbacktrace -o ElmTri

! Output format instructions
200 format (I6, 12E13.3)
300 format (En13.3)
400 format (I6, 10E13.5)
500 format (I7, E19.3)
600 format (I8, I12, I15)
700 format (I10, I16)

120 format (45h Total weight =                             :,F15.3)
122 format (45h Total Sum N*sin(alpha) =                   :,F15.3)
124 format (45h Total Sum N*cos(alpha) =                   :,F15.3)
126 format (45h Total Sum T*sin(alpha) =                   :,F15.3)
128 format (45h Total Sum T*cos(alpha) =                   :,F15.3)
136 format (45h Total Sum U*cos(alpha) =                   :,F15.3)
138 format (45h Total Sum U*sin(alpha) =                   :,F15.3)
140 format (45h Total Sum Vert. loads  =                   :,F15.3)

130 format (45h Sum of All Vert. Forces =                  :,F15.3)
132 format (45h Sum of All Hoz. Forces  =                  :,F15.3)
134 format (55h Sum of all forces and moments imb. in all slices =     :,F5.3)

! Declaration of constants and variables
double precision, parameter :: pi = 4.0*atan(1.0)

double precision Tot_sum_FH, Tot_sum_FV, Tot_sum_N_sin_alpha
double precision Tot_sum_N_cos_alpha, Tot_sum_T_sin_alpha
double precision Tot_sum_T_cos_alpha, old_FS_ave, total_weight
double precision sum_mob_shear, sum_available_shear, FS_ave, factor
double precision slice_function, res_FV, res_FH, res_M
double precision Tot_sum_U_cos_alpha, Tot_sum_U_sin_alpha
double precision Tot_T_res, Tot_T, Tot_vert_loads, min_FS
double precision Tot_sum_C_cos_alpha, Tot_sum_C_sin_alpha
integer i, j, num_slices, num_surfaces, num_forces, num_triangles
integer num_trials, neg_counter, thrust_line_num, critical_surf_num

integer :: output_file=12, input_file=10

double precision, allocatable, dimension (:) :: weight, alpha, beta, Dx
double precision, allocatable, dimension (:) :: slice_alpha, slice_Dx
double precision, allocatable, dimension (:) :: slice_C, slice_phi
double precision, allocatable, dimension (:) :: slice_unit_weight, V_load
double precision, allocatable, dimension (:) :: slice_length, slice_U
double precision, allocatable, dimension (:) :: area, I_vector, L_vector, slice_V_load
double precision, allocatable, dimension (:) :: HL, HR, Ld, N, T, U
double precision, allocatable, dimension (:) :: FS_slice, Total_residue
double precision, allocatable, dimension (:) :: slice_weight, triangle_C
double precision, allocatable, dimension (:) :: triangle_phi, triangle_unit_weight
double precision, allocatable, dimension (:) :: sigma, tau_avail, T_res
double precision, allocatable, dimension (:,:) :: G_matrix
integer, allocatable, dimension (:) :: neg_count_vec

! Input and output file names
open(output_file,file='outputElmTri.txt')
open(input_file,file='inputElmTri.txt')

write(*,*) "--------------------------------------------------------------"
write(*,*) "   ElmTri - Enhanced Limit Method with Triangular slices      "
write(*,*) "--------------------------------------------------------------"
write(*,*) "This program calculates the safety factor of a fully specified"
write(*,*) "   potential slip surface by the enhanced limit method with   "
write(*,*) "          triangular slices by the method proposed by         "
write(*,*) "                 Alexandre and Thomasi (2025)"
write(*,*) "--------------------------------------------------------------"
write(*,*) " VERSION WITH TRIANGULAR END SLICES AND NORMAL FORCE AT DX/3  "
write(*,*) "                                                              "
write(*,*) "    This program is intended for research purposes only.      "
write(*,*) "  No Liability, whether implicit or explicit, is accepted.    "
write(*,*) "--------------------------------------------------------------"
write(*,*) "                                                              "

write(output_file,*) "--------------------------------------------------------------"
write(output_file,*) "   ElmTri - Enhanced Limit Method with Triangular slices      "
write(output_file,*) "--------------------------------------------------------------"
write(output_file,*) "This program calculates the safety factor of a fully specified"
write(output_file,*) "   potential slip surface by the enhanced limit method with   "
write(output_file,*) "          triangular slices by the method proposed by         "
write(output_file,*) "                 Alexandre and Thomasi (2025)"
write(output_file,*) "--------------------------------------------------------------"
write(output_file,*) " VERSION WITH TRIANGULAR END SLICES AND NORMAL FORCE AT DX/3  "
write(output_file,*) "                                                              "
write(output_file,*) "    This program is intended for research purposes only.      "
write(output_file,*) "  No Liability, whether implicit or explicit, is accepted.    "
write(output_file,*) "--------------------------------------------------------------"
write(output_file,*) "                                                              "

read(input_file, *) num_slices

! Calculation of the number of triangles, surfaces and forces
num_triangles = 2 + (num_slices - 2)*2
num_surfaces = 3 + (num_slices - 2)*3
num_forces = 2*num_surfaces

! Number of lines of thrust 
num_trials = 29

write(*,*) "Num slices    Num tria.     Num surfaces"
write(*,600) num_slices, num_triangles, num_surfaces
write(*,*) " "

write(*,*) "Num of unknowns Num of equations"
write(*,700) num_forces, num_triangles*3

write(output_file,*) "Num slices    Num tria.     Num surfaces"
write(output_file,600) num_slices, num_triangles, num_surfaces

write(output_file,*) " "
write(output_file,*) "Num of unknowns Num of equations"
write(output_file,700) num_forces, num_triangles*3
write(output_file,*) " "

! Allocating memory
allocate(weight(num_triangles), area(num_triangles), I_vector(num_forces))
allocate(L_vector(num_forces), HL(num_triangles), HR(num_triangles))
allocate(U(num_triangles), Dx(num_triangles))
allocate(alpha(num_triangles), beta(num_triangles), Ld(num_triangles))
allocate(G_matrix(num_forces, num_forces), N(num_slices))
allocate(T(num_slices), FS_slice(num_slices), neg_count_vec(num_trials))
allocate(Total_residue(num_trials), slice_Dx(num_slices))
allocate(slice_weight(num_slices), slice_alpha(num_slices))
allocate(slice_C(num_slices), slice_phi(num_slices), slice_unit_weight(num_slices))
allocate(triangle_C(num_triangles), triangle_phi(num_triangles))
allocate(triangle_unit_weight(num_triangles))
allocate(sigma(num_slices), tau_avail(num_slices), T_res(num_slices))
allocate(slice_length(num_slices), slice_U(num_slices), V_load(num_triangles))
allocate(slice_V_load(num_slices))

! Initialising vectors 
neg_count_vec = 0
I_vector = 0.0
HL = 0.0
HR = 0.0
Ld = 0.0
U = 0.0

! Initial line of thrust factor as a fraction of the of the "standard"
! line of thrust height (1/3 of the height)
factor = 0.0001

!write(output_file, *) " Initial Factor"
!write(output_file, 300) factor  
    
write(output_file, *) " "

write(output_file, *) " Triangle   Dx          Alpha         Ld          Beta          HL            HR         Area&
                        &          U            c            phi         Gamma        Load"

! Reading Slice width (Dx), slice base angle (alpha) with hor. dir.
! length of the diagonal (Ld), inclined pane (beta), height of the slice
! on the left (HL), height of the slice on the right (HL), area of the
! triangle (Area), pore-pressure resultant at the base of the slice
! cohesion intercept (c), friction angle (phi), unit weight (Gamma)
! and vertical load (Load)
                                          
do i=1, num_triangles
    read(input_file, *) Dx(i), alpha(i), Ld(i), Beta(i), HL(i), HR(i), Area(i), U(i), &
                        & triangle_C(i), triangle_phi(i), triangle_unit_weight(i), V_load(i)
    write(output_file, 200) i, Dx(i), alpha(i), Ld(i), Beta(i), HL(i), HR(i), Area(i), U(i), &
                        & triangle_C(i), triangle_phi(i), triangle_unit_weight(i), V_load(i)
    alpha(i) = alpha(i)*(pi/180)
    beta(i) = beta(i)*(pi/180)
end do

! Total weight of the potential slip surface
total_weight = 0.0

! Calculating weight of each triangle and the total weight
! of the potential slip surface
do i=1, num_triangles
    weight(i) = area(i)*triangle_unit_weight(i)
    total_weight = total_weight + weight(i)
end do

! Assigning soil properties and geometrical dimensions
! to slices
write(*,*) " "
do i=1, num_slices
    
    ! First slice
    if (i==1) then
    
        slice_weight(i) = weight(i)
        
        slice_alpha(i) = alpha(i)
        slice_Dx(i) = Dx(i)
        slice_length(i) = slice_Dx(i)/cos(slice_alpha(i))
        
        slice_U(i) = U(i)
        slice_V_load(i) = V_load(i)
        
        slice_C(i) = triangle_C(i)
        slice_phi(i) = triangle_phi(i)
        slice_unit_weight(i) = triangle_unit_weight(i)
    
    ! Last slice    
    else if (i==num_slices) then
    
        slice_weight(i) = weight(num_triangles)
        
        slice_alpha(i) = alpha(num_triangles)
        slice_Dx(i) = Dx(num_triangles)
        slice_length(i) = slice_Dx(i)/cos(slice_alpha(i))
        
        slice_U(i) = U(num_triangles)
        slice_V_load(i) = V_load(num_triangles)
        
        slice_C(i) = triangle_C(num_triangles)
        slice_phi(i) = triangle_phi(num_triangles)
        slice_unit_weight(i) = triangle_unit_weight(num_triangles)
    
    ! Intermediate slices    
    else
    
        slice_weight(i) = weight(2*i-2) + weight(2*i-1)
        slice_alpha(i) = alpha(2*i-2)
        slice_Dx(i) = Dx(2*i-2)
        slice_length(i) = slice_Dx(i)/cos(slice_alpha(i))
        
        slice_U(i) = U(2*i-2)
        slice_V_load(i) = V_load(2*i-1)
        
        slice_C(i) = triangle_C(2*i-2)
        slice_phi(i) = triangle_phi(2*i-2)
        slice_unit_weight(i) = triangle_unit_weight(2*i-2)
    
    end if

end do

! Variables for identifying the critical line of thrust
thrust_line_num = 0
critical_surf_num = 0
old_FS_ave = 1000.0
min_FS = 1000.0

Total_residue = 0.0

! Loop for the calculation of each line of thrust
do j=1, num_trials
    
    ! Zeroing global matrix and vectors
    
    Tot_T = 0.0                 ! Total tang. force T along slip surface
    Tot_T_res = 0.0             ! Residue of total tang. force T
    Tot_vert_loads = 0.0        ! Total vertical loads applied

    Tot_sum_FH = 0.0            ! Sum of horizontal forces
    Tot_sum_FV = 0.0            ! Sum of vertical forces
    Total_residue = 0.0      ! Total residue force on slip surface

    neg_counter = 0             ! Counter for negative normal forces
    G_matrix = 0.0          ! Global geom. matrix [G]
    L_vector = 0.0                 ! Load vector {L}
    
    sum_mob_shear = 0.0         ! Sum of mobilized shear
    sum_available_shear = 0.0   ! Sum of available shear strength
    
    slice_function = factor     ! Line of thrust multiplier

    ! First slice (first triangle)
    ! Sum of vertical forces
    G_matrix(1,1) = cos(alpha(1))
    G_matrix(1,2) = sin(alpha(1))
    G_matrix(1,4) = -1.0
    L_vector(1) = weight(1) + V_load(1) - U(1)*cos(alpha(1))
    
    ! Sum of horizontal forces
    G_matrix(2,1) = -sin(alpha(1))
    G_matrix(2,2) = cos(alpha(1))
    G_matrix(2,3) = -1
    L_vector(2) = U(1)*sin(alpha(1))
    
    ! Sum of moments
    G_matrix(3,1) = -(Dx(1)/3.0)/cos(alpha(1))
    G_matrix(3,3) = slice_function*Hr(1)/3.0
    L_vector(3) = -weight(1)*Dx(1)/3 + U(1)*(Dx(1)/2.0)/cos(alpha(1)) - V_load(1)*Dx(1)/2.0

    ! Lower triangles
    do i=2, (num_triangles-1), 2

        ! Sum of vertical forces
        G_matrix(3*i-2,3*i-2) = 1.0
        G_matrix(3*i-2,3*i-1) = cos(alpha(i))
        G_matrix(3*i-2,3*i) = sin(alpha(i))
        G_matrix(3*i-2,3*i+1) = -sin(beta(i))
        G_matrix(3*i-2,3*i+2) = -cos(beta(i))
        L_vector(3*i-2) = weight(i) - U(i)*cos(alpha(i))
        
        ! Sum of horizontal forces
        G_matrix(3*i-1,3*i-3) = 1.0
        G_matrix(3*i-1,3*i-1) = -sin(alpha(i))
        G_matrix(3*i-1,3*i) = cos(alpha(i))
        G_matrix(3*i-1,3*i+1) = -cos(beta(i))
        G_matrix(3*i-1,3*i+2) = sin(beta(i))
        L_vector(3*i-1) = U(i)*sin(alpha(i))
        
        ! Sum of moments
        G_matrix(3*i,3*i-2) = -Dx(i)
        G_matrix(3*i,3*i-3) = -(slice_function*HL(i)/3.0 - Dx(i)*tan(alpha(i))) 
        G_matrix(3*i,3*i-1) = -(Dx(i)/cos(alpha(i)))/2.0
        G_matrix(3*i,3*i+1) = slice_function*Ld(i)/3.0
        L_vector(3*i) = -weight(i)*2.0*Dx(i)/3.0 + U(i)*(Dx(i)/2.0)/cos(alpha(i))
        
    
    end do
    
    ! Upper triangles
    do i=3, num_triangles, 2

        !slice_function = factor

        ! Sum of vertical forces
        G_matrix(3*i-2,3*i-2) = sin(beta(i))
        G_matrix(3*i-2,3*i-1) = cos(beta(i))
        G_matrix(3*i-2,3*i+1) = -1
        L_vector(3*i-2) = weight(i) + V_load(i)
        
        ! Sum of horizontal forces
        G_matrix(3*i-1,3*i-2) = cos(beta(i))
        G_matrix(3*i-1,3*i-1) = -sin(beta(i))
        G_matrix(3*i-1,3*i) = -1
        L_vector(3*i-1) = 0
        
        ! Sum of moments
        G_matrix(3*i,3*i-2) = -slice_function*Ld(i)/3.0
        G_matrix(3*i,3*i) = slice_function*HR(i)/3.0
        L_vector(3*i) = -weight(i)*Dx(i)/3.0 - V_load(i)*Dx(i)/2
        
    end do

    ! Last slice (last triangle)
    
    ! Sum of vertical forces
    G_matrix(3*num_triangles-2,3*num_triangles-2) = 1
    G_matrix(3*num_triangles-2,3*num_triangles-1) = cos(alpha(num_triangles))
    G_matrix(3*num_triangles-2,3*num_triangles) = sin(alpha(num_triangles))
    L_vector(3*num_triangles-2) = weight(num_triangles) - U(num_triangles)*cos(alpha(num_triangles)) + V_load(num_triangles)
    
    ! Sum of horizontal forces
    G_matrix(3*num_triangles-1,3*num_triangles-3) = 1
    G_matrix(3*num_triangles-1,3*num_triangles-1) = -sin(alpha(num_triangles))
    G_matrix(3*num_triangles-1,3*num_triangles) = cos(alpha(num_triangles))

    L_vector(3*num_triangles-1) = U(num_triangles)*sin(alpha(num_triangles))
    
    ! Sum of moments
    G_matrix(3*num_triangles,3*num_triangles-3) = -slice_function*HL(num_triangles)/3.0
    G_matrix(3*num_triangles,3*num_triangles-1) = (Dx(num_triangles)/3.0)/cos(alpha(num_triangles))
        
    L_vector(3*num_triangles) = weight(num_triangles)*Dx(num_triangles)/3.0 &
        & - U(num_triangles)*(Dx(num_triangles)/2.0)/cos(alpha(num_triangles)) &
        & + V_load(num_triangles)*(Dx(num_triangles)/2)
    
    ! Solving matrix equation [G]*{i} = {L}
    call gauss_2(G_matrix,L_vector,I_vector,num_forces)
    
    ! Checking equilibrium of first slice
    res_FV = -weight(1) + I_vector(1)*cos(alpha(1)) + I_vector(2)*sin(alpha(1)) - I_vector(4) &
           & + U(1)*cos(alpha(1)) - V_load(1)
           
    res_FH =  -I_vector(1)*sin(alpha(1)) + I_vector(2)*cos(alpha(1)) - I_vector(3) &
           & - U(1)*sin(alpha(1))
           
    res_M = weight(1)*Dx(1)/3.0 - I_vector(1)*(Dx(1)/3.0)/cos(alpha(1)) + &
          & I_vector(3)*slice_function*HR(1)/3.0 - U(1)*(Dx(1)/2)/cos(alpha(1)) &
          &  + V_load(1)*Dx(1)/2.0
          
    Total_residue(j) = Total_residue(j) + abs(res_FV) + abs(res_FH) + abs(res_M)
           
           
    !write(*,*) Total_residue(j)
    !write(*,*) "  "

    ! Checking equilibrium of each triangle
    do i=2, (num_triangles-1), 2
    
        res_FV = -weight(i) + I_vector(3*i-1)*cos(alpha(i)) + I_vector(3*i)*sin(alpha(i)) &
               & + I_vector(3*i-2) -I_vector(3*i+1)*sin(beta(i)) -I_vector(3*i+2)*cos(beta(i)) &
               & + U(i)*cos(alpha(i))
               
        res_FH = I_vector(3*i-3) - I_vector(3*i-1)*sin(alpha(i)) + I_vector(3*i)*cos(alpha(i)) &
               & -I_vector(3*i+1)*cos(beta(i)) + I_vector(3*i+2)*sin(beta(i))  &
               & -U(i)*sin(alpha(i))
               
        res_M = -I_vector(3*i-2)*Dx(i) - I_vector(3*i-3)*(slice_function*HL(i)/3.0 -Dx(i)*tan(alpha(i))) &
               & + weight(i)*Dx(i)*2.0/3.0  - I_vector(3*i-1)*(Dx(i)/2)/cos(alpha(i)) &
               & + I_vector(3*i+1)*slice_function*Ld(i)/3 - U(i)*(Dx(i)/2.0)/cos(alpha(i))
               
                 
        Total_residue(j) = Total_residue(j) + abs(res_FV) + abs(res_FH) + abs(res_M)
        !write(*,*) Total_residue(j)
    
    end do
    !write(*,*) "  "
    
    do i=3, num_triangles, 2
    
        res_FV = -weight(i) + I_vector(3*i-2)*sin(beta(i)) + I_vector(3*i-1)*cos(beta(i)) &
               & -I_vector(3*i+1) - V_load(i)
               
        res_FH = I_vector(3*i-2)*cos(beta(i)) - I_vector(3*i-1)*sin(beta(i)) - I_vector(3*i)
        
        res_M = -I_vector(3*i-2)*slice_function*Ld(i)/3.0 + weight(i)*Dx(i)/3.0 &
              & + I_vector(3*i)*slice_function*HR(i)/3.0 + V_load(i)*Dx(i)/2.0
               
                 
        Total_residue(j) = Total_residue(j) + abs(res_FV) + abs(res_FH) + abs(res_M)   
        !write(*,*) Total_residue(j)
    
    end do
    
    ! check equilibrium of the last slice
    res_FV = -weight(num_triangles) + I_vector(num_forces-2) + &
           & I_vector(num_forces-1)*cos(alpha(num_triangles)) + &
           & I_vector(num_forces)*sin(alpha(num_triangles)) &
           & + U(num_triangles)*cos(alpha(num_triangles)) &
           & - V_load(num_triangles)
           
    res_FH = I_vector(num_forces-3) - I_vector(num_forces-1)*sin(alpha(num_triangles)) &
           & + I_vector(num_forces)*cos(alpha(num_triangles)) &
           & - U(num_triangles)*sin(alpha(num_triangles))
           
    res_M = -weight(num_triangles)*Dx(num_triangles)/3.0 &
          & - I_vector(num_forces-3)*HL(num_triangles)*slice_function/3.0 &
          & + I_vector(num_forces-1)*(Dx(num_triangles)/3.0)/cos(alpha(num_triangles)) &
          & + U(num_triangles)*(Dx(num_triangles)/2)/cos(alpha(num_triangles)) &
          & - V_load(num_triangles)*Dx(num_triangles)/2.0
          
    Total_residue(j) = Total_residue(j) + abs(res_FV) + abs(res_FH) + abs(res_M)
    !write(*,*) "  "
    !write(*,*) Total_residue(j)
    !write(*,*) "  "

    ! Normal force at the base of slice
    ! Tangential force at the base of slice
    do, i=1, num_slices
    
        if (i==1) then
        
            N(1) = I_vector(1)
            T(1) = I_vector(2)
            
        else if (i==num_slices) then
        
            N(num_slices) = I_vector(num_forces-1)
            T(num_slices) = I_vector(num_forces)
        
        else
    
            N(i) = I_vector(6*i-7)     
            T(i) = I_vector(6*i-6)     
            
        end if
    
    end do
    
    ! Counting number of neg. normal forces (N) for each line of thrust
    do i=1, num_slices
    
        if (N(i)<0.0) then
            neg_counter = neg_counter + 1
            
        else 
        
        end if

        sum_mob_shear = sum_mob_shear + T(i)
        sum_available_shear = sum_available_shear + N(i)*tan((pi/180)*slice_phi(i))
        FS_slice(i) = (tan((pi/180)*slice_phi(i)))/(T(i)/N(i))
        
    end do
    
    
    ! Counting number of neg. normal interslice 
    ! forces for each line of thrust
    do i=1, num_forces
    
        if ((mod(i,2)/=0) .and. (I_vector(i)<0.0)) then
            
            neg_counter = neg_counter + 1
            
        else
    
        end if
    
    end do

    ! Total of the negative normal forces for each
    ! line of thrust
    neg_count_vec(j) = neg_counter

    ! Printing line of thrust results
    if (neg_count_vec(j) == 0) then
        
        thrust_line_num = thrust_line_num + 1
        
        write(output_file,*) " "
        write(output_file,*) "------------------------------------"
        write(output_file,*) "Thrust Line Number", thrust_line_num
        write(output_file,*) "------------------------------------"
        
        
        write(output_file,*) " "
        write(output_file,*) "Force number      Force"
        
        N(1) = I_vector(1)
        T(1) = I_vector(2)
        
        do, i=2, num_slices
        
            N(i) = I_vector(6*i-7)
    
        end do
        
        do, i=1, num_forces
        
            write(output_file,500) i, I_vector(i)
    
        end do
        
        write(output_file,*) " "
        write(output_file, *) "  Slice      N            T   &
                            &   Alpha (deg.) Slice length      C      &
                            & Phi (deg.)     Weight      Sigma Eff.  Tau avail.    &
                            & T res.        "
        write(output_file,*)
 
        Tot_sum_FH = 0.0
        
        Tot_sum_N_cos_alpha = 0.0
        Tot_sum_N_sin_alpha = 0.0
        Tot_sum_T_cos_alpha = 0.0
        Tot_sum_T_sin_alpha = 0.0
        
        Tot_sum_U_cos_alpha = 0.0
        Tot_sum_U_sin_alpha = 0.0
        
        Tot_sum_C_cos_alpha = 0.0
        
        Tot_sum_C_sin_alpha = 0.0
        
        
        Tot_T = 0.0
        Tot_T_res = 0.0
        
        Tot_vert_loads = 0.0
        
        do i=1, num_slices
        
            Tot_sum_C_sin_alpha = Tot_sum_C_sin_alpha + &
                                  & slice_C(i)*slice_length(i)*sin(slice_alpha(i))
            Tot_sum_C_cos_alpha = Tot_sum_C_cos_alpha + &
                                  & slice_C(i)*slice_length(i)*cos(slice_alpha(i))
        
            if (i==1) then
                N(i) = I_vector(1)
                T(i) = I_vector(2)
            else
                N(i) = I_vector(6*i-7)
                T(i) = I_vector(6*i-6)
            end if
            
            sigma(i) = N(i)/slice_length(i)
            tau_avail(i) = slice_C(i) + sigma(i)*tan((pi/180)*slice_phi(i))
            T_res(i) = tau_avail(i)*slice_length(i)
            
            N(1) = I_vector(1)
            T(1) = I_vector(2)
            
            sigma(1) = N(1)/slice_length(1)
            tau_avail(1) = slice_C(1) + sigma(1)*tan((pi/180)*slice_phi(i))
            T_res(1) = tau_avail(1)*slice_length(1)
            
            FS_slice(i) = T_res(i)/T(i)
            
            Tot_T = Tot_T + T(i)
            Tot_T_res = Tot_T_res + T_res(i)
            
            Tot_sum_U_cos_alpha = Tot_sum_U_cos_alpha + slice_U(i)*cos(slice_alpha(i))
            Tot_sum_U_sin_alpha = Tot_sum_U_sin_alpha + slice_U(i)*sin(slice_alpha(i))
            
            Tot_sum_N_cos_alpha = Tot_sum_N_cos_alpha + N(i)*cos(slice_alpha(i))
            Tot_sum_N_sin_alpha = Tot_sum_N_sin_alpha + N(i)*sin(slice_alpha(i))
            Tot_sum_T_cos_alpha = Tot_sum_T_cos_alpha + T(i)*cos(slice_alpha(i))
            Tot_sum_T_sin_alpha = Tot_sum_T_sin_alpha + T(i)*sin(slice_alpha(i))
            
            Tot_vert_loads = Tot_vert_loads + slice_V_load(i)

        end do
        
        ! Calculation of the average safety factor
        ! for each line line of thrust
        FS_ave = Tot_T_res/Tot_T
        
        ! Printing results depending whether the 
        ! safety factor is positive or negative
        if (FS_ave > 0.0) then
        
            do i=1, num_slices
            
                write(output_file,400) i, N(i), T(i), slice_alpha(i)*180.0/pi, & 
                                   & slice_length(i), slice_C(i), slice_phi(i), &
                                   & slice_weight(i), sigma(i), tau_avail(i), T_res(i)
        
            end do
 
            if (FS_ave<old_FS_ave) then
                min_FS = FS_ave
                critical_surf_num = thrust_line_num
                old_FS_ave = FS_ave
            else
                min_FS = old_FS_ave
            end if
            
            write(output_file,*) " "
            write(output_file,120) -Total_weight
            write(output_file,124) Tot_sum_N_cos_alpha
            write(output_file,126) Tot_sum_T_sin_alpha
            write(output_file,136) Tot_sum_U_cos_alpha
            write(output_file, 140) -Tot_vert_loads

            write(output_file,*) "-----------------------------------------------------------"
            write(output_file,130) (-Total_weight  - Tot_vert_loads + &
                                    & Tot_sum_N_cos_alpha + Tot_sum_T_sin_alpha &
                                    & + Tot_sum_U_cos_alpha)
            
            write(output_file,*) " "
    
            write(output_file,122) -Tot_sum_N_sin_alpha 
            write(output_file,128) Tot_sum_T_cos_alpha
            write(output_file,138) -Tot_sum_U_sin_alpha
    
            write(output_file,*) "-----------------------------------------------------------"
            write(output_file,132) (-Tot_sum_N_sin_alpha + Tot_sum_T_cos_alpha &
                                    & - Tot_sum_U_sin_alpha)
            
            write(output_file,*) " "
            write(output_file, 134) Total_residue(j)
    
            write(output_file,*) " "
            write(output_file,*) "     Factor"
            write(output_file,300) factor
            
            write(output_file,*) " "
            write(output_file,*) "  Average FS"
            write(output_file,300) FS_ave
            
        else
            write(output_file,*) "Negative safety factor - Discarted line of thrust"
        
        end if
    
    end if
    
    ! incremenenting the line of thrust factor    
    factor = factor + 0.1
	 
end do
write(output_file,*) " "
write(output_file,*) "Thrust Line Number" 
write(output_file,*) critical_surf_num

! Closing input and output files
close(output_file)
close(input_file)

end program ElmTri

!***********************************************************************************************************************************
! SUBROUTINE FOR SOLVING LINEAR SYSTEM
!***********************************************************************************************************************************    
   
  subroutine gauss_2(a,b,x,n)
!===========================================================
! Solutions to a system of linear equations A*x=b
! Method: Gauss elimination (with scaling and pivoting)
! Alex G. (November 2009) for more information please see https://ww2.odu.edu/~agodunov/computing/programs/index.html 
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! b(n)   - array of the right hand coefficients b
! n      - number of equations (size of matrix A)
! output ...
! x(n)   - solutions
! comments ...
! the original arrays a(n,n) and b(n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), b(n), x(n)
double precision s(n)
double precision c, pivot, store
integer i, j, k, l

! step 1: begin forward elimination
do k=1, n-1

! step 2: "scaling"
! s(i) will have the largest element from row i 
  do i=k,n                       ! loop over rows
    s(i) = 0.0
    do j=k,n                    ! loop over elements of row i
      s(i) = max(s(i),abs(a(i,j)))
    end do
  end do

! step 3: "pivoting 1" 
! find a row with the largest pivoting element
  pivot = abs(a(k,k)/s(k))
  l = k
  do j=k+1,n
    if(abs(a(j,k)/s(j)) > pivot) then
      pivot = abs(a(j,k)/s(j))
      l = j
    end if
  end do

! Check if the system has a singular matrix
  if(pivot == 0.0) then
    write(*,*) ' The matrix is singular '
    return
  end if

! step 4: "pivoting 2" interchange rows k and l (if needed)
if (l /= k) then
  do j=k,n
     store = a(k,j)
     a(k,j) = a(l,j)
     a(l,j) = store
  end do
  store = b(k)
  b(k) = b(l)
  b(l) = store
end if

! step 5: the elimination (after scaling and pivoting)
   do i=k+1,n
      c=a(i,k)/a(k,k)
      a(i,k) = 0.0
      b(i)=b(i)- c*b(k)
      do j=k+1,n
         a(i,j) = a(i,j)-c*a(k,j)
      end do
   end do
end do

! step 6: back substitution 
x(n) = b(n)/a(n,n)
do i=n-1,1,-1
   c=0.0
   do j=i+1,n
     c= c + a(i,j)*x(j)
   end do 
   x(i) = (b(i)- c)/a(i,i)
end do

end subroutine gauss_2

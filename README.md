# LES Filtering & Vortical Structures

This program is used for two purpose. 

First it is applied to filter the fully-developed channel DNS data ( Re_tau = 395 ). Gaussian filter is used to filter DNS data.

Second it is used to get vortical structures from three different identification schemes
  - Q criteria ( Haller.1988 )
  - Lambda_2 criteria ( Jeong.1995 )
  - Lambda_ci criteria ( Zhou.1999 )

### 0. Related papers & data 
  - Related paper : https://www.dropbox.com/s/osqskke88qgzow7/Turbulence%20Modeling%20-%20Project2.pdf?dl=0 
  - Related report : https://www.dropbox.com/s/3b180nsfyvw7gkj/report.pdf?dl=0
  - Example data : https://www.dropbox.com/s/241nuubybzb8tgg/instantaneous_velocity_field_re644.plt?dl=0

### 1. Setting for LES_Filtering code
  - Make 'RESULT' folder 
  - Add 3D Velocity data (u,v,w in x,y,z direction) ( tecplot format ) to RESULT folder
  - Change settings in LES_Filtering_setup.f90 code only
    - Set file name variable ( file_name ) 
    - Set the number of meshes ( Nx,Ny,Nz )
    - Set the filter width constant ( FW ) 
    - Set Y-slice point ( N,YP )
    - Set Statistic type ( VS_ONLY )
    - Set vortical structure identification scheme ( VS_CASE )
    - Set the number of filter ( FILTER_OX )
      - If you want to make only vortical structure, set Filter_OX = 0, VS_ONLY = 2
      - If you want to get eddy viscosity and Smargorinsky coefficient, set FILTER_OX =2, VS_ONLY = 0 or 1.
    - Set Y sorting order ( Y_ORDER)
      - Whether the Y is saved in ascending or descending order.
      
    
### 2. Output files
  - Filtered mean velocities ( 3D field, Y slice contour, Y profile )
  - Filtered residual stress ( Y slice contour, Y profile )
  - Filtered strain rate ( Y slice contour, Y profile )
  - Filtered rotation rate ( Y slice contour, Y profile )
  - Filtered eddy viscosity ( Y slice contour, Y profile )
  - Filtered Smargorinsky coefficient ( Y slice contour, Y profile )
  - Vortical Structures 
  
### 3. Code composition
  - LES_filtering_main.f90
  - LES_filtering_module.f90
  - LES_filtering_setup.f90
  - LES_filtering_read.f90
  - LES_filtering_filter.f90
  - LES_filtering_second_filter.f90
  - LES_filtering_output.f90
  - LES_filtering_eig33.f90
  - Vortical_Structures.f90

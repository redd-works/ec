"""
EC9.PY DOCUMENTATION

"""

from inputs import *
from math import sqrt, pi


def sec_class_(beta, f_o=f_o, class_M=class_M):

    """

    :method: material class selection.

    :parameters: beta: width-to-thickness ratio b/t, \n
                f_o: Yield Strength (N/mm^2), \n
                class_M: material of section. \n


    :return: beta: width-to-thickness ratio b/t, \n
             beta1,beta2,beta3 : limits for slenderness parameter, \n
             class_s: Material classification. \n


    :reference: BS EN 1999-1-1 :2007, table 6.2 .           
    """
    # Table 6.2
    epsilon = sqrt(250/f_o)
    beta_1 = epsilon*11
    beta_2 = epsilon*16
    beta_3 = epsilon*22

    class_s = 4
    if class_M == 'A':
        if beta < beta_1:
            class_s = 1
        elif beta < beta_2:
            class_s = 2
        elif beta < beta_3:
            class_s = 3
    return class_s, beta, beta_1, beta_2, beta_3

def sec_class():

    """
    :method: Classifying section based on maximum value returned from 
                using ifferent beta.

    :return: class_s: Material classification. \n


    :reference: BS EN 1999-1-1 :2007, 6.1.4.4 .

    """
    # Section classification
    #### 6.1.4.4
    beta_w = 0.4*b_w/d_w
    class_w, _, _, _, _ = sec_class_(beta_w)

    beta_f = b_f / d_f
    class_f, _, _, _, _ = sec_class_(beta_f)

    class_ = max(class_w, class_f)
    return class_


def comp_buck(f_o=f_o, E=E, A=A, I=I_z, L=L_b, class_M=class_M):

    """

    :method: Calculate design buckling resistance.

    :parameters: f_o: Yield Strength (N/mm^2),\n
                 E: modulus of elasticity (Pa), \n
                 A: Section area (mm^2) \n
                 I_z: Second moment of inertia about z axis (mm^4) \n
                 L: Span(m) \n
                 class_M: class of material from sec_class().\n


    :return: N_rd is the design buckling resistance of the compression member


    :reference: BS EN 1999-1-1 :2007, 6.3.1 .   




    """


    N_pl = f_o*A
    N_cr = pi**2 * E * I / L**2
    lamb_ = sqrt(N_pl/N_cr)
    if class_M == 'A':
        alpha = 0.2
        lamb_0 = 0.1
    phi = 0.5*(1 + alpha*(lamb_ - lamb_0) + lamb_**2)
    xi = 1 / (phi + sqrt(phi**2 - lamb_**2))
    N_rd = xi*N_pl/gamma_m1
    return N_rd


def alpha_3(beta, beta2, beta3, W_el, W_pl):
    return (1 + (beta3-beta)/(beta3-beta2)*(W_pl/W_el - 1))
            

def ltb_buck(f_o=f_o, E=E, C1=1., A=A, h=h, I_y=I_y, I_z=I_z, L=L_d):

    """

    :method: Calculate design buckling resistance moment.

    :parameters: f_o: Yield Strength (N/mm^2),\n
                E: modulus of elasticity (Pa), \n
                A: Section area (mm2) \n
                I_z: Second moment of inertia about z axis (mm^4) \n
                I_y: Second moment of inertia about y axis (mm^4) \n
                h: Depth of section (mm)

    :return: M_rd is the design buckling resistance moment.


    :reference: BS EN 1999-1-1 :2007, 6.3.2.1 .        
    """


    Iw = (I_y * h**2) / 4
    M_cr = C1 * pi**2 * E * I_z / L**2 * (Iw/I_z + L**2*G*J/(pi**2*E*I_z))**0.5
    W_el = I_y/Z_cog
    M_el = f_o*W_el
    W_pl = Wy_plastic()

    beta_w = 0.4*b_w/d_w
    class_, beta, beta1, beta2, beta3 = sec_class_(beta_w)
    if class_ >= 3:
        alpha_3u = alpha_3(beta, beta2, beta3, W_el, W_pl)
        alpha = 0.2
        lamb_0 = 0.4

    lamb_ = sqrt(alpha_3u*M_el/M_cr)
    phi = 0.5*(1 + alpha*(lamb_ - lamb_0) + lamb_**2)
    xi = 1 / (phi + sqrt(phi**2 - lamb_**2))
    M_rd = xi*M_el/gamma_m1

    return M_rd


def Wy_plastic(h=h, Z_cog=Z_cog, t_b=t_f):

    """

    :method: Calculate plastic modulus of section.

    :parameters: Z_cog: centroid of section (mm),\n
                 h: Depth of section (mm),\n
                 t_b = 

    :return: Plastic modulus of section


    :reference:        
    """

    z = h - Z_cog
    h_w = z - t_b
    return 2*(2*(t_w*h_w**2/2) + b*t_b*(z-t_b/2))


def Wz_plastic(A=964, Y_cog=24):

    """

    :method: Calculate plastic modulus of section.

    :parameters: A: area of half-section (mm),
                 Y_cog: centroid of section (mm),\n

    :return: Plastic modulus of section


    :reference:        
    """

    return 2*A*Y_cog


def sec_rd(f_o, E):

    """

    :method: Calculate section resistance, normal, shear, moment and torsional.

    :parameters: f_o: Yield Strength (N/mm^2),\n
                E: modulus of elasticity (Pa), \n
                

    :return: N_rd: normal resistance, \n
                Vy_rd: shear resistance about y axis,\n
                Vz_rd: shear resistance about z axis,\n
                T_rd: Torional resistance, \n
                Myy_rd: moment resistance about y axis, \n
                Mzz_rd: moment resistance about z axis.
         
    """
    class_ = sec_class()
    if class_ == 3:

        # Axial resistance
        N_rd = comp_buck(f_o=f_o, E=E)

        # Moment resistance
        Myy_rd = ltb_buck(f_o=f_o, E=E)
        Wy_pl = Wy_plastic()
        alpha_y = Wy_pl*Z_cog/I_y
        ita_0 = max(min(alpha_y**2, 2), 1)
                       
        Wz_el = I_z/(b/2)
        Wz_pl = Wz_plastic()
        beta_f = b_f/d_f
        class_, beta, beta1, beta2, beta3 = sec_class_(beta_f)
        if class_ >= 2:
            alpha_3u = alpha_3(beta, beta2, beta3, Wz_el, Wz_pl)
        Mzz_rd = alpha_3u*f_o*Wz_el/gamma_m1

        # Shear resistance
        A_v = 2*70*d_f
        Vz_rd = f_o*A_v/(sqrt(3)*gamma_m1)
        A_h = b*2 + b_f*d_f
        Vy_rd = f_o*A_h/(sqrt(3)*gamma_m1)

        # Torsion resistance
        t_max = 6 # mm
        T_rd = f_o*J/(sqrt(3)*gamma_m1*t_max)

                
        return N_rd, Vy_rd, Vz_rd, T_rd, Myy_rd, Mzz_rd, ita_0

    else:
        print('Different class')



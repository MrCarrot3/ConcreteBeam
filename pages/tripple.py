from cmath import sqrt
import math
import streamlit as st
import numpy as np
import plotly.graph_objs as go

#class Interpol:
def linear_LHS( max, len, x):
    m = max/len
    y = m*(x-len)+max
    return y

def mid (lhs, rhs, len, x):
    m=(math.sqrt(lhs*lhs)-math.sqrt(rhs*rhs))/(-len)
    y = m*(x-len)+math.sqrt(rhs*rhs)
    return y

def linear_RHS(max, len, x):
    m = max/len
    y = -m*(x)+max
    return y

def linear (max_neg, len, a, x ):
    b = len - a
    if x<=a:
        m = -max_neg/a
        y = m*x
    else:
        m = (-max_neg-0)/(b)
        y = m*(a-x)-max_neg
    return y


def point_max_mom(l, x, p):
    b = l-x
    moment = p*x*b/l
    return moment

def par_(max, len, x):
    mid = len/2
    a = -max/(mid*mid)
    y = a*(x-mid)*(x-mid)+max
    return y

def udl_max(w, l):
    udl = w*l*l/8
    return udl

def moment_dist(l1, l2 , l3, w1, w2, w3, kn1, lx1, kn2, lx2,kn3, lx3):
    s_1 = 3/l1 #stiffness
    s_2 = 4/l2
    s_3 = 3/l3  # stiffness

    b1=l1-lx1
    b2=l2-lx2
    b3=l3-lx3

    df_a1 = s_1/(s_1+s_2)
    df_a2 = s_2 / (s_1 + s_2)
    df_b1 = s_2 / (s_3 + s_2)
    df_b2 = s_3 / (s_3 + s_2)

    mom1 = w1*l1*l1/8 + kn1*lx1*b1*(l1+b1)/(2*l1*l1)
    mom2_1 = w2*l2*l2/12 + (kn2*lx2*b2*b2)/(l2*l2)
    mom2_2 = w2*l2*l2/12 + (kn2*lx2*lx2*b2)/(l2*l2)
    mom3 = w3*l3*l3/8 + kn3*lx3*b3*(l3+b3)/(2*l3*l3)

                  #0        #1                  #2                          #3                          #4               #5
    r1 = np.array([0,       -mom1,              mom2_1,                   -mom2_2,                          mom3,                   0])
    r2 = np.array([0,       0,                  -(r1[1]+r1[2]),         0,0,0])
    r3 = np.array([0,       r2[2]*df_a1,        r2[2]*df_a2,            (r2[2]*df_a2)/2,                0,                      0])
    r4 = np.array([0,       0,                  0,                      -(r1[3]+r1[4]+r3[3]),           0,                      0])
    r5 = np.array([0,       0,                  0,                      r4[3]*df_b1,                    r4[3]*df_b2,            0])
    r6 = np.array([0,       0,                  r5[3]/2,                0,                              0,                      0])
    r7 = np.array([0,       -r6[2]*df_a1,       -r6[2]*df_a2,           -(r6[2]*df_a2)/2,                0,                     0])
    r8 = np.array([0,       0,                  0,                      -r7[3]*df_b1,                   -r7[3]*df_b2,           0])
    r9 = np.array([0,       r1[1]+r3[1]+r7[1],  r1[2]+r3[2]+r6[2]+r7[2], r1[3]+r3[3]+r5[3]+r7[3]+r8[3], r1[4]+r3[4]+r5[4]+r8[4], 0])
    return r9

def Alpha2 (fc):
    Alpha2_def = 0.85-0.0015*fc
    return Alpha2_def
def Y_conc (fc):
    Y_conc_def = 0.97-0.0025*fc
    return Y_conc_def
def T_Steel (bar_size, no_bars, Fsy):
    Ast = pow(bar_size / 2, 2) * 3.145 * no_bars
    T_kN = Ast*Fsy
    return T_kN
def Moment_Arm (depth, ydn, cover):
    Z_arm = depth-cover-ydn/2
    return Z_arm
def conc_strain (neutral_axis , d_neutral_axis, ec):
    conc_strain = (ec/neutral_axis)*d_neutral_axis
    return conc_strain
def shear_rhs (mom1, mom2, udl, point, l):
    by = -(sqrt(mom1*mom1) - (udl*l*l)/2)/4
    ay = (-udl*l + by)
    r1 = np.array([ay, by])
    return r1
def shear_lhs (mom1, mom2, udl, point, l):       # ^------) moment is on the right hand side
    by = -(sqrt(mom1*mom1) - (udl*l*l)/2)/4
    ay = (-udl*l + by)
    r1 = np.array([ay, by])
    return r1
def shear_rhs (mom1, mom2, udl, point, l ):     # (------^ moment is on the left hand side
    ay = (mom1 + udl*l*(l/2))/4
    by = -(udl*l + ay)
    r1 = np.array([ay, by])
    return r1

def mom_cap(t, z):
    mom = t*z
    return mom

# formatting the page
st.set_page_config(layout="wide")

col1, col2, col3 = st.columns([1,1,2])

with col1:
    User_bar_size = st.selectbox("What type of steel are you using? (mm^2)", (10,12,16,20,24,32,26) )
    User_no_bars = st.slider("How many bars of Steel?", min_value =0, max_value =10, step=1 )
    # use input first, then prompt for the integer using int
    User_fc = st.selectbox("What is the grade of Concrete? (MPa)", (20, 25, 32, 40, 50 ))
    User_Width = st.slider("What is the width of Concrete? (mm)", min_value=0, max_value=2000, step = 5)
    User_Depth = st.slider("What is the depth of Concrete? (mm)", min_value=0, max_value=2000, step = 5)
    User_cover = st.slider("What is the cover of Concrete? (mm)", min_value= 0 , max_value= 100, step=5)

#print(T_Steel(User_bar_size, User_no_bars, 500)) #this needs the steel inputs which are AST & Fsy
# neutral axis calc
tension_in_steel = T_Steel(User_bar_size, User_no_bars, 500)
dn = tension_in_steel/(Alpha2(User_fc)*User_fc*Y_conc(User_fc)*User_Width)
ydn = dn*Y_conc(User_fc)
strain_steel = conc_strain(dn,User_Depth-dn-User_cover, 0.003)
moment_capacity = mom_cap(tension_in_steel/1000,
                           Moment_Arm(User_Depth, ydn, User_cover)/1000)

#print(round(dn, 2),"mm is the depth of your Neutral axis")
with col2:
    len1 = st.slider("length of 1st span", min_value=1, max_value=10, step=1 )
    len2 = st.slider("length of 2nd Span", min_value=1, max_value=10, step=1 )
    len3 = st.slider("length of 3rd Span", min_value=1, max_value=10, step=1 )
    udl1 = st.slider("Select the Uniformly distributed load for span 1 ", min_value=0, max_value=20, step=1 )
    udl2 = st.slider("Select the Uniformly distributed load for span 2", min_value=0, max_value=20, step=1 )
    udl3 = st.slider("Select the Uniformly distributed load for span 3", min_value=0, max_value=20, step=1 )

with (col3):
    kn1= st.slider("Select the point load for span 1", min_value=0, max_value=100, step=10 )
    lx1 = st.slider("Select the x coordinate  for span 1", min_value=1, max_value=len1, step=1 )
    kn2= st.slider("Select the point load for span 2", min_value=0, max_value=100, step=10 )
    lx2 = st.slider("Select the x coordinate  for span 2", min_value=1, max_value=len2, step=1 )
    kn3= st.slider("Select the point load for span 3", min_value=0, max_value=100, step=10 )
    lx3 = st.slider("Select the x coordinate  for span 3", min_value=1, max_value=len3, step=1 )

area = 3.145*User_bar_size*User_bar_size/4
effective_width = User_Width - User_cover

if User_bar_size == 1:
    barx = np.array([User_Width/2])
else:
    barx = np.linspace(User_cover, effective_width, User_no_bars)
    bary = np.full_like(barx, User_cover)

beam_display = go.Figure()
beam_display.add_trace(go.Scatter(x=[0, 0, User_Width, User_Width], y=[0, User_Depth,User_Depth, 0], fill = "toself"))
beam_display.add_trace(go.Scatter(
    x = barx,
    y = bary,
    mode = "markers",
    marker = dict(
        size = User_bar_size,
        sizemode = "diameter",
    )
    )
)
beam_display.update_layout(showlegend=False)
beam_display.update_xaxes(range=[0, 2200])
beam_display.update_yaxes(range=[0, 2200])

st.plotly_chart(beam_display)

print(round(Moment_Arm(User_Depth,ydn, User_cover), 2),"kNm")

print()

if strain_steel > 0.0025:
    print(f"Steel has strained at {strain_steel: 3f}")
else:
    print(f"Steel not has strained at {strain_steel: 3f}")

#######Main MOMENT DIST######
max_moment = moment_dist(len1, len2, len3, udl1, udl2, udl3, kn1, lx1, kn2, lx2,kn3, lx3)
st.write(max_moment)

len1_max_point = point_max_mom(len1, lx1, kn1)
len2_max_point = point_max_mom(len2, lx2, kn2)
len3_max_point = point_max_mom(len3, lx3, kn3)

mom_x = np.array([0,
            len1/10,
            2*len1/10,
            3*len1/10,
            4*len1/10,
            5*len1/10,
            6*len1/10,
            7*len1/10,
            8*len1/10,
            9*len1/10,
            len1,
            len1+len2 / 10,
            len1+2 * len2 / 10,
            len1+3 * len2 / 10,
            len1+4 * len2 / 10,
            len1+5 * len2 / 10,
            len1+6 * len2 / 10,
            len1+7 * len2 / 10,
            len1+8 * len2 / 10,
            len1+9 * len2 / 10,
            len1+len2,
            len1+len2+len3/10,
            len1+len2+2*len3/10,
            len1+len2+3*len3/10,
            len1+len2+4*len3/10,
            len1+len2+5*len3/10,
            len1+len2+6*len3/10,
            len1+len2+7*len3/10,
            len1+len2+8*len3/10,
            len1+len2+9*len3/10,
            len1 + len2 +len3])
mom_y_udl = np.array([0,
        -linear_LHS(max_moment[1], len1, len1/10 )-               par_(udl_max(udl1, len1), len1, len1/10),
        -linear_LHS(max_moment[1], len1, 2*len1/10) -             par_(udl_max(udl1, len1), len1, 2*len1/10),
        -linear_LHS(max_moment[1], len1, 3*len1/10)-              par_(udl_max(udl1, len1), len1, 3*len1/10),
        -linear_LHS(max_moment[1], len1, 4 * len1 / 10) -         par_(udl_max(udl1, len1), len1, 4 * len1 / 10),
        -linear_LHS(max_moment[1], len1, 5 * len1 / 10) -         par_(udl_max(udl1, len1), len1, 5 * len1 / 10),
        -linear_LHS(max_moment[1], len1, 6*len1 / 10) -           par_(udl_max(udl1, len1), len1, 6*len1 / 10),
        -linear_LHS(max_moment[1], len1, 7 * len1 / 10)-          par_(udl_max(udl1, len1), len1, 7 * len1 / 10),
        -linear_LHS(max_moment[1], len1, 8 * len1 / 10) -         par_(udl_max(udl1, len1), len1, 8 * len1 / 10),
       - linear_LHS(max_moment[1], len1, 9 * len1 / 10) -         par_(udl_max(udl1, len1), len1, 9 * len1 / 10),
        -max_moment[1],

         #middle
         mid(max_moment[2],max_moment[3], len2, len2/10)-          par_(udl_max(udl2, len2), len2, len2/10),
         mid(max_moment[2],max_moment[3], len2, 2*len2 / 10) -     par_(udl_max(udl2, len2), len2, 2*len2 / 10),
         mid(max_moment[2],max_moment[3], len2, 3 * len2 / 10)-    par_(udl_max(udl2, len2), len2, 3 * len2 / 10),
         mid(max_moment[2],max_moment[3], len2, 4 * len2 / 10) -   par_(udl_max(udl2, len2), len2, 4 * len2 / 10),
         mid(max_moment[2],max_moment[3], len2, 5 * len2 / 10) -   par_(udl_max(udl2, len2), len2, 5 * len2 / 10),
         mid(max_moment[2],max_moment[3], len2, 6 * len2 / 10) -   par_(udl_max(udl2, len2), len2, 6 * len2 / 10),
         mid(max_moment[2],max_moment[3], len2, 7 * len2 / 10) -   par_(udl_max(udl2, len2), len2, 7 * len2 / 10),
         mid(max_moment[2],max_moment[3], len2, 8 * len2 / 10)-    par_(udl_max(udl2, len2), len2, 8 * len2 / 10),
         mid(max_moment[2],max_moment[3], len2, 9 * len2 / 10) -   par_(udl_max(udl2, len2), len2, 9 * len2 / 10),
         -max_moment[3],

        #right hand side
         linear_RHS(max_moment[4], len3, len3/10)-           par_(udl_max(udl3, len3), len3, len3/10),
         linear_RHS(max_moment[4], len3, 2*len3 / 10) -      par_(udl_max(udl3, len3), len3, 2*len3 / 10),
         linear_RHS(max_moment[4], len3, 3 * len3 / 10)-     par_(udl_max(udl3, len3), len3, 3 * len3 / 10),
         linear_RHS(max_moment[4], len3, 4 * len3 / 10) -    par_(udl_max(udl3, len3), len3, 4 * len3 / 10),
         linear_RHS(max_moment[4], len3, 5 * len3 / 10) -    par_(udl_max(udl3, len3), len3, 5 * len3 / 10),
         linear_RHS(max_moment[4], len3, 6 * len3 / 10) -    par_(udl_max(udl3, len3), len3, 6 * len3 / 10),
         linear_RHS(max_moment[4], len3, 7 * len3 / 10) -    par_(udl_max(udl3, len3), len3, 7 * len3 / 10),
         linear_RHS(max_moment[4], len3, 8 * len3 / 10)-     par_(udl_max(udl3, len3), len3, 8 * len3 / 10),
         linear_RHS(max_moment[4], len3, 9 * len3 / 10) -    par_(udl_max(udl3, len3), len3, 9 * len3 / 10),
         linear_RHS(max_moment[4], len3, 10 * len3 / 10)-    par_(udl_max(udl3, len3), len3, 10 * len3 / 10),

         ])
mom_y_point = np.array([0,
        linear(len1_max_point, len1, lx1, 1*len1/10),
        linear(len1_max_point, len1, lx1, 2*len1/10),
        linear(len1_max_point, len1, lx1, 3*len1/10),
        linear(len1_max_point, len1, lx1, 4*len1/10),
        linear(len1_max_point, len1, lx1, 5*len1/10),
        linear(len1_max_point, len1, lx1, 6*len1/10),
        linear(len1_max_point, len1, lx1, 7*len1/10),
        linear(len1_max_point, len1, lx1, 8*len1/10),
        linear(len1_max_point, len1, lx1, 9*len1/10),
        linear(len1_max_point, len1, lx1, 10*len1/10),
        #middel
        linear(len2_max_point, len2, lx2, 1*len2/10),
        linear(len2_max_point, len2, lx2, 2*len2/10),
        linear(len2_max_point, len2, lx2, 3*len2/10),
        linear(len2_max_point, len2, lx2, 4*len2/10),
        linear(len2_max_point, len2, lx2, 5*len2/10),
        linear(len2_max_point, len2, lx2, 6*len2/10),
        linear(len2_max_point, len2, lx2, 7*len2/10),
        linear(len2_max_point, len2, lx2, 8*len2/10),
        linear(len2_max_point, len2, lx2, 9*len2/10),
        linear(len2_max_point, len2, lx2, 10*len2/10),
        #right hand side
        linear(len3_max_point, len3, lx3, 1*len3/10),
        linear(len3_max_point, len3, lx3, 2*len3/10),
        linear(len3_max_point, len3, lx3, 3*len3/10),
        linear(len3_max_point, len3, lx3, 4*len3/10),
        linear(len3_max_point, len3, lx3, 5*len3/10),
        linear(len3_max_point, len3, lx3, 6*len3/10),
        linear(len3_max_point, len3, lx3, 7*len3/10),
        linear(len3_max_point, len3, lx3, 8*len3/10),
        linear(len3_max_point, len3, lx3, 9*len3/10),
        linear(len3_max_point, len3, lx3, 10*len3/10),
         ])
mom_tot = mom_y_point + mom_y_udl

moment_plot = go.Figure()
moment_plot.add_trace(go.Scatter(
    x = mom_x,
    y = mom_tot))


st.write(mid(max_moment[2],max_moment[3], len2, 5*len2/10)- par_(udl_max(udl1, len2), len2, 5*len2/10))

moment_plot.add_trace(go.Scatter(
    x = (0, len2+len1+len3),
    y = (moment_capacity, moment_capacity)

    ))
st.plotly_chart(moment_plot)






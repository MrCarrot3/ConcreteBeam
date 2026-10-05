from cmath import sqrt
#from main import poi_mom_rhs

import math
import streamlit as st
import numpy as np
import plotly.express as px
import plotly.graph_objs as go

#class Interpol:
def linear_LHS( max, len, x):
    m = max/len
    y = m*(x-len)+max
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
def moment_dist(l1, l2 , w1, w2, a1, a2, f1, f2):
    s_1 = 3/l1 #stiffness
    s_2 = 3/l2

    b1 = l1 - a1
    b2 = l2 - a2

    df_1 = s_1/(s_1+s_2)
    df_2 = s_2/(s_1 + s_2)

    mom1 = (w1*l1*l1)/8 + f1*a1*b1*(l1+b1)/(2*l1*l1)  # add point load moment
    mom2 = (w2*l2*l2)/8 + f2*a2*b2*(l2+b2)/(2*l2*l2)  # add point load moment

    r1 = np.array([0, -mom1, mom2, 0])
    r2 = np.array([0, (r1[1]+r1[2])*df_1, (r1[1]+r1[2])*df_2, 0])
    r3 = np.array([0, -(-r1[1]+r2[1]), r1[2]-r2[2], 0])

    return r3
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
def spans_x(start, l, lx, n = 100 ):
    local = np.unique(np.concatenate([np.linspace(start, l, n+1),[lx+start]]))
    #linespace will return the values from start to finsih,
    # however your actual point load locaiton of x is not in this list,
    # therefore you will need to add it
    #this is where concatinate comes in, it will add the locaiton of point load
    #lx
    return local
def spans_y(start, l, lx, n = 100 ):
    local = np.unique(
            np.concatenate(
                [np.linspace(start, l, n+1),[lx+start]]
            )
    )
    #linespace will return the values from start to finsih,
    # however your actual point load locaiton of x is not in this list,
    # therefore you will need to add it
    #this is where concatinate comes in, it will add the locaiton of point load
    #lx
    localremoved = local[:-1]
    return localremoved
#def mom_cap(t, z):
    #mom = t*z
    #return mom
def poi_mom_rhs(l, f, a):
    b = l-a
    r1 = np.array([f*a*b*(l+b)/(2*l*l),
                   (f*a*a*b*(2*l+b)/(2*l*l*l))]
                  )
    return r1
def double_mom(b, d, cover_t, cover_, bars_t, size_t, bars_c, size_c, alpha, y_conc, fc):
    # tensile
    ats = 3.145 * (size_t / 2) * (size_t / 2) * bars_t
    print(f"Area of tensile Steel {ats} mm2")
    ts = ats * 500 / 1000
    print(f"Tensile strength {ts} kN")

    # compressive steel
    asc = 3.145 * (size_c / 2) * (size_c / 2) * bars_c
    print(f"Area of compressive steel {asc} mm2")

    a = y_conc*(fc/1000)*alpha*b
    b_ = ((200 * 1000) * asc * 0.003 / 1000 - ts)
    c = -((200*1000)*asc* 0.003/1000)*cover_t
    print(f"{a} dn  {b_}dn + {c}dn")

    dn = ((-b_)+ math.sqrt(b_*b_ - 4*a*c))/(2*a)
    print(f"The neutral axis is located at: {dn}mm")
    print(f"Area of steel{ats} mm2")
    print(f"Tensile strength{ts} kN")

    #mu = ts1*(d-dcs) + ts2*(d - y*dn/2)
    # ts1 = esc * E * asc
    esc = 0.003*(dn-cover_t)/dn
    #print(f"cover bottom is{cover_b} mm")
    ts1 = esc*(200*1000)*(asc/1000)
    #print(f"strength of compressive steel in {ts1} kN")
    ts2 = a*dn
    print(f"Concrete Compressive strennght is {ts2}kN")

    mu = (ts1*(d-cover_t) + ts2*(d - y_conc*dn/2))/1000
    return mu
def mom_cap(b, d, cover_b, bars_t, size_t, fc ):
    t = (size_t/2)*(size_t/2)*3.145*bars_t*500/1000
    print(t)
    alpha = 0.85 - 0.0015*fc
    y = 0.97 - 0.0025*fc
    dn  = t/(y*alpha*fc*b/1000)
    print(dn)
    z = d - y*dn/2 - cover_b
    mom = t*z/1000
    return mom


# formatting the page
st.set_page_config(layout="wide")

col1, col2, col3 = st.columns([1,1,2])

with col1:
    st.header("Input Concrete section elements")
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
#moment_capacity = mom_cap(tension_in_steel/1000,
                           #Moment_Arm(User_Depth, ydn, User_cover)/1000)

#print(round(dn, 2),"mm is the depth of your Neutral axis")
with col2:
    st.header("Input Spans & Loads")
    len1 = st.slider("length of 1st span", min_value=0, max_value=10, step=1 )
    len2 = st.slider("length of 2nd Span", min_value=0, max_value=10, step=1 )
    udl1 = st.slider("Select the Uniformly distributed load  for Span 1 ", min_value=0, max_value=20, step=1 )
    udl2 = st.slider("Select the Uniformly distributed load for Span 2", min_value=0, max_value=20, step=1 )
    kn1= st.slider("Select the point load for span 1", min_value=0, max_value=100, step=1 )
    lx1 = st.slider("Select the x coordinate  for span 1", min_value=1, max_value=len1, step=1 )
    kn2= st.slider("Select the point load for span 2", min_value=0, max_value=100, step=1 )
    lx2 = st.slider("Select the x coordinate for span 2", min_value=1, max_value=len2, step=1 )

with (col3):
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
#print(round(Moment_Arm(User_Depth,ydn, User_cover), 2),"kNm")

#print()

if strain_steel > 0.0025:
    print(f"Steel has strained at {strain_steel: 3f}")
else:
    print(f"Steel not has strained at {strain_steel: 3f}")

#THIS IS THE DEFINING MOMENT OF THE  DOCUMENT
max_moment = moment_dist(len1, len2, udl1, udl2, lx1 , lx2, kn1, kn2)
len1_max_point = point_max_mom(len1, lx1, kn1)
len2_max_point = point_max_mom(len2, lx2, kn2)


#st.write(max_moment)
#max_moment_1 = max(max_moment[1], max_moment[2])
#max_moment_2 = -max_moment[2]
max_moment_1 = udl_max(udl1, len1)
max_moment_2 = udl_max(udl2, len2)


#setting up the graph

x1 = spans_x(0, len1, lx1) #x for span1 values distributed throughout
x2 = spans_x(len1, len1+ len2, lx2) #x for span1 + span2 values distributed throughout


mom_x = np.unique((np.concatenate([x1,x2])))
len1_mom = spans_y(0, len1, lx1) #x for span1 values for moment calc
len2_mom = spans_x(0, len2, lx2) #x for span3 values for moment calc
                                        # its weird, but for the last span
                                        #instead of using span_y, use
                                        # span_x instead because it works better
                                        # it honestly for got the reason
                                        # but the graph breaks if you don't
y_spans = (np.concatenate([len1_mom, len2_mom]))

mom_udl_len1 = [] #<---- this is the array that the moment for length 1 goes into
mom_udl_len2 = [] #<---- this is the array that the moment for length 2 goes into

q = 0
for q in range(len(len1_mom)):
    mom1_udl = -linear_LHS(max_moment[1], len1, len1_mom[q]) - par_(udl_max(udl1, len1), len1, len1_mom[q]) + linear(len1_max_point, len1, lx1, len1_mom[q])
    mom_udl_len1.append(float(mom1_udl)) #<---- append the moment from each point into the array for length 1
    q = q+1 #<-- -increase the increment by 1
#st.write(mom_udl_len1) <---- FOR DEBUGGING

#span 2 moment

j = 0
for j in range(len(len2_mom)):
    mom3_udl = (linear_RHS(max_moment[2], len2, len2_mom[j]) - par_(udl_max(udl2, len2), len2, len2_mom[j]) + linear(len2_max_point, len2, lx2, len2_mom[j])
                )
    mom_udl_len2.append(float(mom3_udl))
    j = j+1

mom_y_udl = np.concatenate([mom_udl_len1, mom_udl_len2[1:]])




# shear

rhs_shear = poi_mom_rhs(len1, kn1,  lx1) #
st.write(f" The shear forces on the right hand side is  {rhs_shear}")

#mom_x = np.array([0,
#            len1/10,
#            2*len1/10,
#            3*len1/10,
#            4*len1/10,
#            5*len1/10,
#            6*len1/10,
#            7*len1/10,
#            8*len1/10,
#            9*len1/10,
#            len1,
#            len1+len2 / 10,
#            len1+2 * len2 / 10,
#            len1+3 * len2 / 10,
#            len1+4 * len2 / 10,
#            len1+5 * len2 / 10,
#            len1+6 * len2 / 10,
#            len1+7 * len2 / 10,
#            len1+8 * len2 / 10,
#            len1+9 * len2 / 10,
#            len1+len2])

#mom_y_udl = np.array([0,
#        linear_LHS(max_moment[1], len1, len1/10 )-
#         par_(udl_max(udl1, len1), len1, len1/10),
#        linear_LHS(max_moment[1], len1, 2*len1/10) -
#         par_(udl_max(udl1, len1), len1, 2*len1/10),
#        linear_LHS(max_moment[1], len1, 3*len1/10)-
#         par_(udl_max(udl1, len1), len1, 3*len1/10),
#         linear_LHS(max_moment[1], len1, 4 * len1 / 10) -
#         par_(udl_max(udl1, len1), len1, 4 * len1 / 10),
#         linear_LHS(max_moment[1], len1, 5 * len1 / 10) -
#         par_(udl_max(udl1, len1), len1, 5 * len1 / 10),
#         linear_LHS(max_moment[1], len1, 6*len1 / 10) -
#         par_(udl_max(udl1, len1), len1, 6*len1 / 10),
#         linear_LHS(max_moment[1], len1, 7 * len1 / 10)-
#         par_(udl_max(udl1, len1), len1, 7 * len1 / 10),
#         linear_LHS(max_moment[1], len1, 8 * len1 / 10) -
#         par_(udl_max(udl1, len1), len1, 8 * len1 / 10),
#         linear_LHS(max_moment[1], len1, 9 * len1 / 10) -
#         par_(udl_max(udl1, len1), len1, 9 * len1 / 10),
#         max(max_moment[2], max_moment[1]),
#         #right hand side
#         linear_RHS(max_moment[2], len2, len2/10)- par_(udl_max(udl2, len2), len2, len2/10),
#         linear_RHS(max_moment[2], len2, 2*len2 / 10) - par_(udl_max(udl2, len2), len2, 2*len2 / 10),
#         linear_RHS(max_moment[2], len2, 3 * len2 / 10)-par_(udl_max(udl2, len2), len2, 3 * len2 / 10),
#         linear_RHS(max_moment[2], len2, 4 * len2 / 10) - par_(udl_max(udl2, len2), len2, 4 * len2 / 10),
#         linear_RHS(max_moment[2], len2, 5 * len2 / 10) - par_(udl_max(udl2, len2), len2, 5 * len2 / 10),
#         linear_RHS(max_moment[2], len2, 6 * len2 / 10) - par_(udl_max(udl2, len2), len2, 6 * len2 / 10),
#         linear_RHS(max_moment[2], len2, 7 * len2 / 10) - par_(udl_max(udl2, len2), len2, 7 * len2 / 10),
#         linear_RHS(max_moment[2], len2, 8 * len2 / 10)- par_(udl_max(udl2, len2), len2, 8 * len2 / 10),
#         linear_RHS(max_moment[2], len2, 9 * len2 / 10) - par_(udl_max(udl2, len2), len2, 9 * len2 / 10),
#         linear_RHS(max_moment[2], len2, 10 * len2 / 10)- par_(udl_max(udl2, len2), len2, 10 * len2 / 10),
#         ])
#mom_y_point = np.array([0,
#        linear(len1_max_point, len1, point_loc_1, 1*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 2*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 3*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 4*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 5*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 6*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 7*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 8*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 9*len1/10),
#        linear(len1_max_point, len1, point_loc_1, 10*len1/10),
#        #right hand side
#        linear(len2_max_point, len2, point_loc_2, 1*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 2*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 3*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 4*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 5*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 6*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 7*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 8*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 9*len2/10),
#        linear(len2_max_point, len2, point_loc_2, 10*len2/10),
#         ])

#mom_tot = mom_y_point + mom_y_udl

mom_tot = mom_y_udl


moment_plot = go.Figure()
moment_plot.add_trace(go.Scatter(
    x = mom_x,
    y = mom_tot))


#moment_plot.add_trace(go.Scatter(
#    x = (0, len2+len1),
#    y = (moment_capacity, moment_capacity)

#    ))

#####Double moment adding####
no_sections= int(st.number_input("How many sections do you want?", step= 1))
moment_cap_x=[]
moment_cap_y=[]
for i in range(0,no_sections):
    col1, col2 = st.columns(2)
    with col1:

        form = st.number_input("From", key=f"in1{i}", step=0.5)
        moment_cap_x.append(form)
        to = st.number_input("To", key=f"in2{i}", step=0.5)
        moment_cap_x.append(to)
    with col2:
        ck = st.checkbox("Steel at top", key=f"in8{i}")
        if ck == True:
            ck = st.checkbox("Compression Steel?", key=f"in9{i}")
            if ck == True:
                moment = -double_mom(User_Width,
                                    User_Depth,
                                    User_cover,
                                    User_cover,
                                    st.number_input("Number of bottom bars", key=f"in7{i}", step=1),
                                    st.selectbox("Size of bottom bars(mm^2)",
                                                                (10, 12, 16, 20, 24, 32, 26), key=f"in4{i}"),
                                    st.number_input("Number of Top bars", key=f"in5{i}", step=1),
                                    st.selectbox("Size of top bars(mm^2)",
                                                (10, 12, 16, 20, 24, 32, 26), key=f"in6{i}"),
                                     0.8083,
                                    0.85,
                                    User_fc)
            else:
                moment = - mom_cap(User_Width,
                                   User_Depth,
                                   User_cover,
                                   st.number_input("Number of bottom bars", key=f"in10{i}", step=1),
                                   st.selectbox("Size of bottom bars(mm^2)",
                                                (10, 12, 16, 20, 24, 32, 26), key=f"in11{i}"),
                                   User_fc)
        else:
            ck = st.checkbox("Compression Steel?", key=f"in9{i}")
            if ck == True:
                moment = double_mom(User_Width,
                                User_Depth,
                                User_cover,
                                User_cover,
                                st.number_input("Number of bottom bars", key=f"in7{i}", step=1),
                                st.selectbox("Size of bottom bars(mm^2)",
                                            (10, 12, 16, 20, 24, 32, 26), key=f"in4{i}"),
                                st.number_input("Number of Top bars", key=f"in5{i}", step=1),
                                st.selectbox("Size of top bars(mm^2)",
                                            (10, 12, 16, 20, 24, 32, 26), key=f"in6{i}"),                                    0.8083,
                                0.85,
                                User_fc )
            else:
                moment = mom_cap(User_Width,
                                User_Depth,
                                User_cover,
                                st.number_input("Number of bottom bars", key=f"in10{i}", step=1),
                                 st.selectbox("Size of bottom bars(mm^2)",
                                              (10, 12, 16, 20, 24, 32, 26), key=f"in11{i}"),                                User_fc)

        moment_cap_y.append(moment)
        moment_cap_y.append(moment)

moment_plot.add_trace(go.Scatter(
    x =moment_cap_x,
    y =moment_cap_y, )) ####adding the double cross-section moment


moment_plot.add_trace(go.Scatter(
    x = mom_x,
    y = mom_tot))

st.plotly_chart(moment_plot)
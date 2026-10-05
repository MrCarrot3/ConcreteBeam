from cmath import sqrt
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
def mid (lhs, rhs, len, x):
    m=(lhs-rhs)/(0-len)
    y = -m*(len-x)+rhs
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
        m = (-max_neg-0)/b
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

                  #0        #1                  #2                          #3                              #4                              #5
    r1 = np.array([0,       -mom1,              mom2_1,                   -mom2_2,                          mom3,                           0])
    r2 = np.array([0,       0,                  -(r1[1]+r1[2]),             0,                              0,                              0])
    r3 = np.array([0,       r2[2]*df_a1,        r2[2]*df_a2,                (r2[2]*df_a2)/2,                0,                              0])
    r4 = np.array([0,       0,                  0,                          -(r1[3]+r1[4]+r3[3]),           0,                              0])
    r5 = np.array([0,       0,                  0,                          r4[3]*df_b1,                    r4[3]*df_b2,                    0])
    r6 = np.array([0,       0,                  r5[3]/2,                    0,                              0,                              0])
    r7 = np.array([0,       -r6[2]*df_a1,       -r6[2]*df_a2,               -(r6[2]*df_a2)/2,                0,                             0])
    r8 = np.array([0,       0,                  0,                          -r7[3]*df_b1,                   -r7[3]*df_b2,                   0])
    r9 = np.array([0,       (r1[1]+r3[1]+r7[1]),  (r1[2]+r3[2]+r6[2]+r7[2]), (r1[3]+r3[3]+r5[3]+r7[3]+r8[3]), (r1[4]+r3[4]+r5[4]+r8[4]),      0])
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

col1, col2, col3, col4 = st.columns([1,1,1,2])

with col1:
    st.header("Input Concrete section elements")
    User_bar_size = st.selectbox("What type of steel are you using? (mm^2)", (10,12,16,20,24,32,26) )
    User_no_bars = st.slider("How many bars of Steel?", min_value =0, max_value =10, step=1 )
    # use input first, then prompt for the integer using int
    User_fc = st.selectbox("What is the grade of Concrete? (MPa)", (20, 25, 32, 40, 50 ))
    User_Width = st.slider("What is the width of Concrete? (mm)", min_value=0, max_value=1000, step = 10)
    User_Depth = st.slider("What is the depth of Concrete? (mm)", min_value=0, max_value=1000, step = 10)
    User_cover = st.slider("What is the cover of Concrete? (mm)", min_value= 0 , max_value= 100, step=5)

with col2:
    st.header("Input Span lengths & UDLs")
    len1 = int(st.number_input("length of 1st span"))
    len2 = int(st.number_input("length of 2nd Span",  step=1 ))
    len3 = int(st.number_input("length of 3rd Span",  step=1 ))
    udl1 = int(st.number_input("Select the Uniformly distributed load for span 1 ", min_value=0, step=1 ))
    udl2 = int(st.number_input("Select the Uniformly distributed load for span 2", min_value=0, step=1 ))
    udl3 = int(st.number_input("Select the Uniformly distributed load for span 3", min_value=0, step=1 ))

with (col3):
    st.header("Input Point loads & locations")
    kn1= int(st.number_input("Select the point load for span 1", min_value=0, max_value=100, step=10 ))
    lx1 = int(st.number_input("Select the x coordinate  for span 1", min_value=1, max_value=len1, step=1 ))
    kn2= int(st.number_input("Select the point load for span 2", min_value=0, max_value=100, step=10 ))
    lx2 = int(st.number_input("Select the x coordinate  for span 2", min_value=1, max_value=len2, step=1 ))
    kn3= int(st.number_input("Select the point load for span 3", min_value=0, max_value=100, step=10 ))
    lx3 = int(st.number_input("Select the x coordinate  for span 3", min_value=1, max_value=len3, step=1 ))

with (col4):
    st.header("Concrete beam section")
    area = 3.145 * User_bar_size * User_bar_size / 4
    effective_width = User_Width - User_cover

    if User_bar_size == 1:
        barx = np.array([User_Width / 2])
    else:
        barx = np.linspace(User_cover, effective_width, User_no_bars)
        bary = np.full_like(barx, User_cover)

    beam_display = go.Figure()
    beam_display.add_trace(
        go.Scatter(x=[0, 0, User_Width, User_Width], y=[0, User_Depth, User_Depth, 0], fill="toself"))
    beam_display.add_trace(go.Scatter(
        x=barx,
        y=bary,
        mode="markers",
        marker=dict(
            size=User_bar_size,
            sizemode="diameter",
            color = 'red'
            )
    )
    )
    beam_display.update_layout(showlegend=False)
    beam_display.update_xaxes(range=[0, 1000])
    beam_display.update_yaxes(range=[0, 2000])

    st.plotly_chart(beam_display)

#######Main MOMENT DIST######
max_moment = moment_dist(len1, len2, len3, udl1, udl2, udl3, kn1, lx1, kn2, lx2,kn3, lx3)
#   st.write(max_moment)

#st.write(f"Moment dist is {max_moment}")  <---- FOR DEBUGGING
len1_max_point = point_max_mom(len1, lx1, kn1) #max point load at span 1
len2_max_point = point_max_mom(len2, lx2, kn2) #max point load at span 2
len3_max_point = point_max_mom(len3, lx3, kn3) #max point load at span 3

#setting up the graph

x1 = spans_x(0, len1, lx1) #x for span1 values distributed throughout
x2 = spans_x(len1, len1+ len2, lx2) #x for span1 + span2 values distributed throughout
x3 = spans_x(len2+len1, len1+len2+len3, lx3) #x for span1 + span2+  span3 values distributed throughout

mom_x = np.unique((np.concatenate([x1,x2,x3])))
#st.write(mom_x) <---- FOR DEBUGGING

#st.write(f"elgnths for span 1 {x1}")   <---- FOR DEBUGGING
#st.write(f"elgnths for span 2 {x2}")   <---- FOR DEBUGGING
#st.write(f"elgnths for span 3 {x3}")   <---- FOR DEBUGGING

# ****** THIS IS SETTING UP ALL THE MOMENT PLOTS FOR INDIVIDUAL SPANS *****
# THE X VALUES FOR INDIVIDUAL SPAN CALCS HAVE TO MATCH THE X VALUES FOR THE TOTAL X SPAN
# E.G.      x PLOT          0 - 1 - 2 -3 -4 -5 -6 - 7 -8 -9 -10
#      X PLOT FOR Y CALC    0 - 1 - 2 -3 -0 -1 -2  -3 -4 -0 -1
len1_mom = spans_y(0, len1, lx1) #x for span1 values for moment calc
len2_mom = spans_y(0, len2, lx2) #x for span2 values for moment calc
len3_mom = spans_x(0, len3, lx3) #x for span3 values for moment calc
                                        # its weird, but for the last span
                                        # instead of using span_y, use
                                        # span_x instead because it works better
                                        # it honestly for got the reason
                                        # but the graph breaks if you don't
y_spans = (np.concatenate([len1_mom, len2_mom, len3_mom]))


mom_udl_len1 = [] #<---- this is the array that the moment for length 1 goes into
mom_udl_len2 = [] #<---- this is the array that the moment for length 2 goes into
mom_udl_len3 = [] #<---- this is the array that the moment for length 3 goes into

#st.write(max_moment) <---- more debugging


########## FOR DEBUGGING ################ This tells you if the number of x for the plot values
#if len(mom_x) == len(y_spans):         # equals the number of x values for the moment calc
#    st.write("All Goods")              # which will calculate the y values.
#    st.write(f"x {len(mom_x)} Spans")  # You need the correct x values to calculate the y values
#    st.write(f"y {len(y_spans)} Spans")#
#else:                                  #
#    st.write("No Bueno")               #
#    st.write(f"x {len(mom_x)} Spans")  #
#    st.write(f"y {len(y_spans)} Spans")#
#########################################

#span 1 moment
q = 0
for q in range(len(len1_mom)):
    mom1_udl = -linear_LHS(max_moment[1], len1, len1_mom[q]) - par_(udl_max(udl1, len1), len1, len1_mom[q]) + linear(len1_max_point, len1, lx1, len1_mom[q])
    mom_udl_len1.append(float(mom1_udl)) #<---- append the moment from each point into the array for length 1
    q = q+1 #<-- -increase the increment by 1
#st.write(mom_udl_len1) <---- FOR DEBUGGING

#span 2 moment
w = 0
for w in range(len(len2_mom)):
    mom2_udl = (mid(max_moment[2], max_moment[4], len2, len2_mom[w]) - par_(udl_max(udl2, len2), len2, len2_mom[w]) + linear(len2_max_point, len2, lx2, len2_mom[w])
                )
    mom_udl_len2.append(float(mom2_udl))
    w = w+1
#st.write(mom_udl_len2) <---- FOR DEBUGGING

#span 3 moment
j = 0
for j in range(len(len3_mom)):
    mom3_udl = (linear_RHS(max_moment[4], len3, len3_mom[j]) - par_(udl_max(udl3, len3), len3, len3_mom[j]) + linear(len3_max_point, len3, lx3, len3_mom[j])
                )
    mom_udl_len3.append(float(mom3_udl))
    j = j+1
#st.write(mom_udl_len3) <---- FOR DEBUGGING
#st.write(f"Number of elements in the y array is {len(len1_mom), len(len2_mom), len(len3_mom)}"
         #f" Number of elements in the x arrary is {len(x1), len(x2), len(x3)}") <---- FOR DEBUGGING

mom_y_udl = np.concatenate([mom_udl_len1, mom_udl_len2[1:], mom_udl_len3[1:]])

mom_tot = mom_y_udl

moment_plot = go.Figure()


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
                moment = double_mom(User_Width,
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
                moment = mom_cap(User_Width,
                                   User_Depth,
                                   User_cover,
                                   st.number_input("Number of bottom bars", key=f"in10{i}", step=1),
                                   st.selectbox("Size of bottom bars(mm^2)",
                                                (10, 12, 16, 20, 24, 32, 26), key=f"in11{i}"),
                                   User_fc)
        else:
            ck = st.checkbox("Compression Steel?", key=f"in9{i}")
            if ck == True:
                moment = - double_mom(User_Width,
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
                moment = - mom_cap(User_Width,
                                User_Depth,
                                User_cover,
                                st.number_input("Number of bottom bars", key=f"in10{i}", step=1),
                                 st.selectbox("Size of bottom bars(mm^2)",
                                              (10, 12, 16, 20, 24, 32, 26), key=f"in11{i}"),                                User_fc)

        moment_cap_y.append(moment)
        moment_cap_y.append(moment)

##### ADDING IN THE MOMENT CAPACITY GRAPH
moment_plot.add_trace(go.Scatter(
    x =moment_cap_x,
    y =moment_cap_y,
    name = "Moment Capacity (kNm)",)) ####adding the double cross-section moment


#### ADDING IN THE MOMENT FROM THE LOADS
moment_plot.add_trace(go.Scatter(
    x = mom_x,
    y = mom_tot))

### drawing in the locaiton of the supports
moment_plot.add_trace(go.Scatter(
    x = [0,len1, len1+len2, len1+len2+len3],
    y = [0, 0,0,0],

))

moment_plot.update_layout(
    xaxis_title = "(m)", yaxis_title = "Moment (kNm)",
)
#st.write(mid(max_moment[2],max_moment[3], len2, 5*len2/10)- par_(udl_max(udl1, len2), len2, 5*len2/10))
###^^^^ IDK why i wrote this, not sure what I was debugging
st.plotly_chart(moment_plot)

#moment_plot.add_trace(go.Scatter(
#    x = (0, len2+len1+len3),
#    y = (moment_capacity, moment_capacity)

#    ))

def shear_left(m1, w, l, p, px):
    vb1 = m1-w*l*l/2-p*px
    va1 = w*l+p+vb1
    r1 = np.array([va1, vb1])
    return r1


def shear_mid(m1, m2, w, l, p, px):
    va1 = (m1+w*l*l/2 -p*px-m2)/l
    vb1 = (-w*l+va1-p)
    r1 = np.array([va1, vb1])
    return

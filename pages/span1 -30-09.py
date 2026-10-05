import streamlit as st
import numpy as np
import plotly.express as px
import plotly.graph_objs as go

def moment (l, w, p, loc,x ):
    b = l - loc

    if x <= loc:
        mom = w*x*(l-x)/2 + p*b*x/l
        return mom
    else:
        maxmoment = p*loc*b/l
        grad_ = -maxmoment/b
        mom = w*x*(l-x)/2 + (maxmoment+grad_*(x-loc))
        return mom

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

def equivl (j, low, l):
    wp = 8*j*low/(l*l)
    return wp

def mom_cap(t, z):
    mom = t*z
    return mom
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





# formatting the page
st.set_page_config(layout="wide")

col1, col2, col3 = st.columns([1,1,4])

with col1:
    User_bar_size = st.selectbox("What type of steel are you using? (mm^2)", (10,12,16,20,24,32,26) )
    User_no_bars = st.slider("How many bars of Steel?", min_value =0, max_value =10, step=1 )
    # use input first, then prompt for the integer using int
    User_fc = st.selectbox("What is the grade of Concrete? (MPa)", (20, 25, 32, 40, 50 ))
    User_Width = st.slider("What is the width of Concrete? (mm)", min_value=0, max_value=2000, step=10)
    User_Depth = st.slider("What is the depth of Concrete? (mm)", min_value=0, max_value=2000, step=10)
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
    len1 = int(st.number_input("length of beam", min_value=1, max_value=100, step=1 ))
    udl1 = int(st.number_input("Select the Uniformly distributed load ", min_value=0, max_value=1000, step=1 ))
    kn1 = int(st.number_input("Select the point load", min_value=0, max_value=1000, step=1 ))
    lx1 = int(st.number_input("Select the x coordinate", min_value=0, max_value=len1, step=1 ))

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

print(round(Moment_Arm(User_Depth,ydn, User_cover), 2),"kNm")

print()

if strain_steel > 0.0025:
    print(f"Steel has strained at {strain_steel: 3f}")
else:
    print(f"Steel not has strained at {strain_steel: 3f}")

st.write(f"The Tension of Steel is  {T_Steel(User_bar_size, User_no_bars, 500)}")
st.write(f"The Neutral axis is located at {dn}")
moment_plot = go.Figure()

len1_max_point = point_max_mom(len1, lx1, kn1) #max point load at span 1

x1 = spans_x(0, len1, lx1) #x for span1 values distributed throughout
len1_mom = spans_y(0, len1, lx1) #x for span1 values for moment calc


mom_udl_len1 = [] #<---- this is the array that the moment for length 1 goes into

q = 0
for q in range(len(len1_mom)):
    mom1_udl = moment(len1, udl1, kn1, lx1, len1_mom[q])
    mom_udl_len1.append(float(mom1_udl)) #<---- append the moment from each point into the array for length 1
    q = q+1 #<-- -increase the increment by 1
#st.write(mom_udl_len1) <---- FOR DEBUGGING






moment_plot.add_trace(go.Scatter(
    x = x1,
    y =mom_udl_len1
    ),
    )
moment_plot.add_trace(go.Scatter(
    x = (0, len1),
    y = (moment_capacity, moment_capacity),

    ))
st.plotly_chart(moment_plot)


shear_plot = go.Figure()
shear_plot.add_trace(go.Scatter(
    x = (0,
            len1/10,
            2*len1/10,
            3*len1/10,
            4*len1/10,
            5*len1/10,
            6*len1/10,
            7*len1/10,
            8*len1/10,
            9*len1/10,
            len1),
    y = (

    )
))



import matplotlib.pyplot as plt

d_climb = 157.14
d_descent = 158.88
d_cruise = 1850 - d_climb - d_descent 

h_cruise0 = 8000
h_cruise1 = 9000
y_flight_profile = [0, h_cruise0, h_cruise1, 0]
x_flight_profile = [0, d_climb, d_climb+d_cruise, d_climb+d_cruise+d_descent]
y_cruise_step = [h_cruise0, h_cruise0, h_cruise1, h_cruise1]
x_cruise_step = [d_climb, d_climb+d_cruise/2, d_climb+d_cruise/2, d_climb+d_cruise]

plt.plot(x_flight_profile, y_flight_profile)
plt.plot(x_cruise_step, y_cruise_step, linestyle=':', color='black')
plt.xlim([0, 1850])
plt.ylim([0, 10000])
plt.xlabel("Range [km]", size='large')
plt.ylabel("Altitude [m]", size='large')
plt.show()
from src.frame2d.frame2d import *

# Luo tyhjä kehä
frame = Frame2D(num_elements=4)
# Luo pilari (koordinaatit, profile=vapaaehtoinen)
col = SteelColumn([[0,0], [0, 5000]], profile='Ipe 300')
# Lisää pilari kehälle
frame.add(col)
# Lisää jäykkä tuki pisteeseen (0,0)
frame.add(FixedSupport([0, 0]))
# Lisää pistekuorma pilarin yläpäähän
frame.add(PointLoad([0, 5000], [10e3, -200e3, 0]))
# Luo kehän fem -malli
frame.generate()
# Laske
frame.calculate()
# Piirtää kehän (pilarin) näytölle
frame.plot()

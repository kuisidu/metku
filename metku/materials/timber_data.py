"""
Timber material data

Author: Viktor Haimi
"""

from enum import Enum
try:
    from framefem.framefem import Material
except:
    from metku.framefem.framefem import Material


class T(Enum):
    # Timber types
    C14 = 1
    C16 = 2
    C18 = 3
    C20 = 4
    C22 = 5
    C24 = 6
    C27 = 7
    C30 = 8
    C35 = 9
    C40 = 10
    D30 = 11
    D35 = 12
    D40 = 13
    D50 = 14
    D60 = 15
    D70 = 16
    GL24c = 17
    GL24h = 18
    GL28c = 19
    GL28h = 20
    GL30c = 21
    GL30h = 22
    GL32c = 23
    GL32h = 24
    Kerto_S = 25
    Kerto_T = 26
    Kerto_Q_21_24 = 27
    Kerto_Q_27_69 = 28

# solid_timber = sahatavara
# lvl = LVL
# glt = Liimapuu


solid_timber = {T.C14: {"fmk": 14, "ft0k": 8, "ft90k": 0.4, "fc0k": 16, "fc90k": 2.0, "fvk": 1.7, "E0mean": 7000,
                        "E005": 4700, "E90mean": 230, "Gmean": 440, "G005": 300, "rhok": 290, "rhomean": 350},
                T.C16: {"fmk": 16, "ft0k": 10, "ft90k": 0.4, "fc0k": 17, "fc90k": 2.2, "fvk": 2.0, "E0mean": 8000,
                        "E005": 5400, "E90mean": 270, "Gmean": 500, "G005": 340, "rhok": 310, "rhomean": 370},
                T.C18: {"fmk": 18, "ft0k": 11, "ft90k": 0.5, "fc0k": 18, "fc90k": 2.2, "fvk": 2.0, "E0mean": 9000,
                        "E005": 6000, "E90mean": 300, "Gmean": 560, "G005": 380, "rhok": 320, "rhomean": 380},
                T.C20: {"fmk": 20, "ft0k": 12, "ft90k": 0.5, "fc0k": 19, "fc90k": 2.3, "fvk": 2.3, "E0mean": 9500,
                        "E005": 6400, "E90mean": 320, "Gmean": 590, "G005": 405, "rhok": 330, "rhomean": 390},
                T.C22: {"fmk": 22, "ft0k": 13, "ft90k": 0.5, "fc0k": 20, "fc90k": 2.4, "fvk": 2.4, "E0mean": 10000,
                        "E005": 6700, "E90mean": 330, "Gmean": 630, "G005": 430, "rhok": 340, "rhomean": 410},
                T.C24: {"fmk": 24, "ft0k": 14, "ft90k": 0.5, "fc0k": 21, "fc90k": 2.5, "fvk": 2.5, "E0mean": 11000,
                        "E005": 7400, "E90mean": 370, "Gmean": 690, "G005": 460, "rhok": 350, "rhomean": 420},
                T.C27: {"fmk": 27, "ft0k": 16, "ft90k": 0.5, "fc0k": 22, "fc90k": 2.6, "fvk": 2.7, "E0mean": 11500,
                        "E005": 7700, "E90mean": 380, "Gmean": 720, "G005": 480, "rhok": 370, "rhomean": 450},
                T.C30: {"fmk": 30, "ft0k": 18, "ft90k": 0.6, "fc0k": 23, "fc90k": 2.7, "fvk": 3.0, "E0mean": 12000,
                        "E005": 8000, "E90mean": 400, "Gmean": 750, "G005": 500, "rhok": 380, "rhomean": 460},
                T.C35: {"fmk": 35, "ft0k": 21, "ft90k": 0.6, "fc0k": 25, "fc90k": 2.8, "fvk": 3.4, "E0mean": 13000,
                        "E005": 8700, "E90mean": 430, "Gmean": 810, "G005": 540, "rhok": 400, "rhomean": 480},
                T.C40: {"fmk": 40, "ft0k": 24, "ft90k": 0.6, "fc0k": 25, "fc90k": 2.8, "fvk": 3.8, "E0mean": 14000,
                        "E005": 9400, "E90mean": 470, "Gmean": 880, "G005": 590, "rhok": 420, "rhomean": 500},
                T.D30: {"fmk": 30, "ft0k": 18, "ft90k": 0.6, "fc0k": 23, "fc90k": 8.0, "fvk": 4, "E0mean": 11000,
                         "E005": 9200, "E90mean": 730, "Gmean": 690, "G005": 500, "rhok": 530, "rhomean": 640},
                T.D35: {"fmk": 35, "ft0k": 21, "ft90k": 0.6, "fc0k": 25, "fc90k": 8.1, "fvk": 4, "E0mean": 12000,
                         "E005": 10100, "E90mean": 800, "Gmean": 750, "G005": 550, "rhok": 540, "rhomean": 650},
                T.D40: {"fmk": 40, "ft0k": 24, "ft90k": 0.6, "fc0k": 26, "fc90k": 8.3, "fvk": 4, "E0mean": 13000,
                         "E005": 10900, "E90mean": 860, "Gmean": 810, "G005": 600, "rhok": 550, "rhomean": 660},
                T.D50: {"fmk": 50, "ft0k": 30, "ft90k": 0.6, "fc0k": 29, "fc90k": 9.3, "fvk": 4, "E0mean": 14000,
                        "E005": 11800, "E90mean": 930, "Gmean": 880, "G005": 650, "rhok": 620, "rhomean": 750},
                T.D60: {"fmk": 60, "ft0k": 36, "ft90k": 0.6, "fc0k": 32, "fc90k": 10.5, "fvk": 4.5, "E0mean": 17000,
                        "E005": 14300, "E90mean": 1130, "Gmean": 1060, "G005": 700, "rhok": 700, "rhomean": 840},
                T.D70: {"fmk": 70, "ft0k": 42, "ft90k": 0.6, "fc0k": 34, "fc90k": 13.5, "fvk": 5, "E0mean": 20000,
                        "E005": 16800, "E90mean": 1330, "Gmean": 1250, "G005": 750, "rhok": 900, "rhomean": 1080}
                }

glt = {T.GL24c: {"fmk": 24, "ft0k": 14, "ft90k": 0.35, "fc0k": 21, "fc90k": 2.4, "fvk": 2.2, "E0mean": 11600,
                 "E005": 9400, "E90mean": 320, "Gmean": 590, "G005": 480, "rhok": 350, "rhomean": 390},
       T.GL24h: {"fmk": 24, "ft0k": 16.5, "ft90k": 0.4, "fc0k": 24, "fc90k": 2.7, "fvk": 2.7, "E0mean": 11600,
                 "E005": 9400, "E90mean": 390, "Gmean": 720, "G005": 530, "rhok": 380, "rhomean": 390},
       T.GL28c: {"fmk": 28, "ft0k": 16.5, "ft90k": 0.4, "fc0k": 24, "fc90k": 2.7, "fvk": 2.7, "E0mean": 12600,
                 "E005": 10200, "E90mean": 390, "Gmean": 720, "G005": 580, "rhok": 380, "rhomean": 430},
       T.GL28h: {"fmk": 28, "ft0k": 19.5, "ft90k": 0.45, "fc0k": 26.5, "fc90k": 3.0, "fvk": 3.2, "E0mean": 12600,
                 "E005": 10200, "E90mean": 420, "Gmean": 780, "G005": 630, "rhok": 410, "rhomean": 460},
       T.GL30c: {"fmk": 30, "ft0k": 20, "ft90k": 0.5, "fc0k": 25, "fc90k": 3.0, "fvk": 3.5, "E0mean": 13000,
                 "E005": 10800, "E90mean": 300, "Gmean": 650, "G005": 540, "rhok": 390, "rhomean": 430},
       T.GL30h: {"fmk": 30, "ft0k": 19.5, "ft90k": 0.5, "fc0k": 24.5, "fc90k": 2.5, "fvk": 3.5, "E0mean": 13000,
                 "E005": 10800, "E90mean": 300, "Gmean": 650, "G005": 540, "rhok": 390, "rhomean": 430},
       T.GL32c: {"fmk": 32, "ft0k": 19.5, "ft90k": 0.45, "fc0k": 26.5, "fc90k": 3.0, "fvk": 3.2, "E0mean": 13700,
                 "E005": 11100, "E90mean": 420, "Gmean": 780, "G005": 630, "rhok": 410, "rhomean": 470},
       T.GL32h: {"fmk": 32, "ft0k": 22.5, "ft90k": 0.5, "fc0k": 29, "fc90k": 3.3, "fvk": 3.8, "E0mean": 13700,
                 "E005": 11100, "E90mean": 460, "Gmean": 850, "G005": 690, "rhok": 430, "rhomean": 500}
       }

lvl = {T.Kerto_S: {"fmk": 44, "s": 0.12, "fm0flatk": 50, "ft0k": 35, "ft90edgek": 0.8, "fc0k": 35, "fc90edgek": 6.0,
                   "fc90flatk": 1.8, "fvk": 4.1, "fr0k": 2.3, "E0mean": 13800, "E005": 11600, "Gedgemean": 600,
                   "Gedge005": 400, "rhok": 480, "rhomean": 510},
       T.Kerto_T: {"fmk": 27, "s": 0.15, "fm0flatk": 32, "ft0k": 24, "ft90edgek": 0.5, "fc0k": 26, "fc90edgek": 4.0,
                   "fc90flatk": 1.0, "fvk": 2.4, "fr0k": 1.3, "E0mean": 10000, "E005": 8800, "Gedgemean": 400,
                   "Gedge005": 300, "rhok": 410, "rhomean": 440},
       T.Kerto_Q_21_24: {"fmk": 28, "s": 0.12, "fm0flatk": 32, "ft0k": 19, "ft90edgek": 6.0, "fc0k": 19, "fc90edgek": 9.0,
                   "fc90flatk": 2.2, "fvk": 4.5, "fr0k": 1.3, "E0mean": 10000, "E005": 8300, "Gedgemean": 600,
                   "Gedge005": 400, "rhok": 480, "rhomean": 510},
       T.Kerto_Q_27_69: {"fmk": 32, "s": 0.12, "fm0flatk": 36, "ft0k": 26, "ft90edgek": 6.0, "fc0k": 26, "fc90edgek": 9.0,
                   "fc90flatk": 2.2, "fvk": 4.5, "fr0k": 1.3, "E0mean": 10500, "E005": 8800, "Gedgemean": 600,
                   "Gedge005": 400, "rhok": 480, "rhomean": 510}
       }


class Timber(Material):
    def __init__(self, timber_type: T=T.C14, hardness: str='softwood'):
        self.timber_type = timber_type
        self.hardness = hardness
        if timber_type in solid_timber:
            self.type = 'solid_timber'
            self.fmk = solid_timber[timber_type]["fmk"]
            self.ft0k = solid_timber[timber_type]["ft0k"]
            self.ft90k = solid_timber[timber_type]["ft90k"]
            self.fc0k = solid_timber[timber_type]["fc0k"]
            self.fc90k = solid_timber[timber_type]["fc90k"]
            self.fvk = solid_timber[timber_type]["fvk"]
            self.E0mean = solid_timber[timber_type]["E0mean"]
            self.E005 = solid_timber[timber_type]["E005"]
            self.E90mean = solid_timber[timber_type]["E90mean"]
            self.Gmean = solid_timber[timber_type]["Gmean"]
            self.G005 = solid_timber[timber_type]["G005"]
            self.rhok = solid_timber[timber_type]["rhok"]
            self.rhomean = solid_timber[timber_type]["rhomean"]

        elif timber_type in glt:
            self.type = 'glt'
            self.fmk = glt[timber_type]["fmk"]
            self.ft0k = glt[timber_type]["ft0k"]
            self.ft90k = glt[timber_type]["ft90k"]
            self.fc0k = glt[timber_type]["fc0k"]
            self.fc90k = glt[timber_type]["fc90k"]
            self.fvk = glt[timber_type]["fvk"]
            self.E0mean = glt[timber_type]["E0mean"]
            self.E005 = glt[timber_type]["E005"]
            self.E90mean = glt[timber_type]["E90mean"]
            self.Gmean = glt[timber_type]["Gmean"]
            self.G005 = glt[timber_type]["G005"]
            self.rhok = glt[timber_type]["rhok"]
            self.rhomean = glt[timber_type]["rhomean"]

        elif timber_type in lvl:
            self.type = 'lvl'
            self.direction = 'edge'
            self.fmk = lvl[timber_type]["fmk"]
            self.s = lvl[timber_type]["s"]
            self.fm0flatk = lvl[timber_type]["fm0flatk"]
            self.ft0k = lvl[timber_type]["ft0k"]
            self.ft90edgek = lvl[timber_type]["ft90edgek"]
            self.fc0k = lvl[timber_type]["fc0k"]
            self.fc90edgek = lvl[timber_type]["fc90edgek"]
            self.fc90flatk = lvl[timber_type]["fc90flatk"]
            self.fvk = lvl[timber_type]["fvk"]
            self.fr0k = lvl[timber_type]["fr0k"]
            self.E0mean = lvl[timber_type]["E0mean"]
            self.E005 = lvl[timber_type]["E005"]
            self.Gedgemean = lvl[timber_type]["Gedgemean"]
            self.Gmean = self.Gedgemean
            self.Gedge005 = lvl[timber_type]["Gedge005"]
            self.G005 = self.Gedge005
            self.rhok = lvl[timber_type]["rhok"]
            self.rhomean = lvl[timber_type]["rhomean"]

        if timber_type in lvl:
            super().__init__(self.E0mean, self.rhomean, G=self.Gedgemean)
        else:
            super().__init__(self.E0mean, self.rhomean, G=self.Gmean)

    def __repr__(self):
        return f"{str(self.timber_type)[2:]}"

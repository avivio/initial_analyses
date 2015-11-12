x <- read.table(text = '
2.20E-28
1.56E-16
                2.57E-33
                2.57E-33
                1.45E-35
                2.07E-23
                7.10E-35
                4.45E-49
                3.07E-13
                5.27E-03
                2.37E-01
                8.74E-02
                1.00E+00
                4.50E-03
                1.34E-21
                8.96E-01
                1.68E-15
                2.57E-02
                6.94E-09
                1
                ')
x <- as.numeric(x[,1])
p.adjust(x, "holm")
writeClipboard(as.character(p.adjust(x, "holm")))
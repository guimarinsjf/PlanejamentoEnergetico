from cvxopt.modeling import op, variable, matrix
import cvxopt.solvers as solvers

class otimiza(object):
    n_ut = 2
    n_uh = 1
    cdef = 500
    carga = 50
    produtib = 0.95
    vmin = 20
    vmax =  100
    custo_termica = [ 10, 25]
    afluencias = [ [ 23, 16], [19, 14], [15, 11]]
    estagios = 3
    gmax_termica = [ 15, 10]
    engolimento = 60

    def __init__(self):
        self.engolimento = self.engolimento

    def pde(self, discret):

        intervalo = ((self.vmax - self.vmin))/(discret-1)

        VI = []
        volume = self.vmax
        for i in range(discret):
            VI.append(volume)
            volume = volume - intervalo

        # Define variáveis de decisão
        vf = variable(1, 'vf')
        vt = variable(1, 'vt')
        vv = variable(1, 'vv')
        gt = variable(2, 'gt')
        deficit = variable(1, 'deficit')
        alpha =variable(1, 'alpha')

        # Define função objetivo
        objetivo = 0
        objetivo = objetivo +  0.001*vv[0]
        objetivo += self.custo_termica[0]*gt[0]
        objetivo += self.custo_termica[1]*gt[1]
        objetivo += self.cdef*deficit[0]
        objetivo += alpha[0]

        # Define Restricoes
        restricoes = []
        # Atendimento a Demanda
        restricoes.append(self.produtib*vt[0] + gt[0] + gt[1] + deficit[0] == self.carga)
        # Limites ou Canalizacao das Variavies
        restricoes.append(vf[0] >= self.vmin)
        restricoes.append(vf[0] <= self.vmax)
        restricoes.append(vt[0] >= 0)
        restricoes.append(vt[0] <= self.engolimento)
        restricoes.append(vv[0] >= 0)
        restricoes.append(gt[0] >= 0)
        restricoes.append(gt[0] <= self.gmax_termica[0])
        restricoes.append(gt[1] >= 0)
        restricoes.append(gt[1] <= self.gmax_termica[1])
        restricoes.append(deficit[0] >= 0)
        restricoes.append(alpha[0] >= 0)

        media = []
        lagrange = []
        beta =[]
        for estagio in range(self.estagios-1,-1,-1):
            media.insert(0,[])
            lagrange.insert(0,[])
            beta.insert(0,[])
            for idiscret in range(discret):
                media[0].append(0)
                lagrange[0].append(0)
                beta[0].append(0)
                for icen in range(2):

                    # Balanco Hidrico
                    restricoes.append(vf[0] == VI[idiscret] +
                                      self.afluencias[estagio][icen] -
                                      vt[0] - vv[0])
                                      
                                      
                                      
#                    if estagio != self.estagios-1:
#                        for rest in range(len(lagrange[1])):
#                            
#                            
#                            restricoes.append(alpha[0]-lagrange[1][rest]*vf->=0)
                                      
                                      

                    problema = op(objetivo, restricoes)

                    problema.solve('dense', 'glpk')
                    
                   

                    media[0][-1] = media[0][-1] + objetivo.value()[0]
                    lagrange[0][-1] += restricoes[-1].multiplier.value[0]
                    del restricoes[-1]
                    
                    
                media[0][-1] = media[0][-1]/2
                lagrange[0][-1] = lagrange[0][-1]/2
                beta[0][-1] = media[0][-1] - lagrange[0][-1]*VI[idiscret]  
            vasco = 1000


        vasco = 1000

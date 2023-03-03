import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier

from utils import *

class HypoClassifier:
    def __init__(self, hypo, kernel, C, cv=5):
        assert C > 0.0
        self.hypo = hypo
        self.cells_A = hypo.nonclonal_cells
        self.cells_B = hypo.clonal_cells
        assert not (set(self.cells_A) & set(self.cells_B))
        
        self.cv_buckets = cv          # number of cross-correlation parts
        self.kernel = kernel
        self.C = C            # SVM regularization parameter
        self.X, self.y = self.prepare_Xy(self.hypo.batch.quant)
        
        # classifier info
        self.clf, self.title, self.scores = None, None, None
        if min(len(self.cells_A), len(self.cells_B)) <= MIN_CELLS_IN_CLONE:
            print(red('Warning:'),
                  'No learning is done! Classes for {} too small: {} and {} cells'.format(
                  blue(self.hypo.clone_tag()), green(len(self.cells_A)), green(len(self.cells_B))))
            return
        
        self.svc = self.get_svc(kernel, C)
        if self.X.shape[0] > 0 and self.X.shape[1] > 0:
            self.clf = self.learn()
            if cv > 0:
                self.scores = self.dishonest_cross_validate()
        
    def __repr__(self):
        res = ""
        res += "{} features, kernel={}, C={};".format(
               blue(self.X.shape[1]), blue(self.kernel), blue('%.2e' % self.C))
        if not self.scores is None:
            res += (red(" Dishonest acc: ") + green("%0.2f (+/-%0.2f)" %
                   (self.scores.mean(), self.scores.std() * 2)))
        return res

    def prepare_Xy(self, quant):
        X_df = quant.subge(self.cells_A + self.cells_B).ge
        self.X = np.array(X_df)
        self.y = [0] * len(self.cells_A) + [1] * len(self.cells_B)
        return (self.X, self.y)

    def get_svc(self, kernel, C):
        if kernel == 'linear':
            svc = svm.SVC(kernel='linear', C=C, probability=True)
            title = 'SVC with linear kernel'
        elif kernel == 'rbf':
            svc = svm.SVC(kernel='rbf', gamma=0.7, C=C)
            title = 'SVC with RBF kernel'
        elif kernel == 'poly':
            svc = svm.SVC(kernel='poly', degree=3, C=C)
            title = 'SVC with polynomial (degree 3) kernel'
        elif kernel == 'lin':
            svc = svm.LinearSVC(C=C)
            title = 'LinearSVC (linear kernel)'
        elif kernel == 'log-linear':
            svc = LogisticRegression(C=C, penalty='l1', tol=0.01)
            title = 'Logistic Regression'
        elif kernel == 'decision-tree':
            svc = DecisionTreeClassifier(random_state=0)
            title = 'Decision Tree'
        else:
            assert False, 'no such kernel: ' + kernel
        return svc
    
    def learn(self):
        #print("learn over {} x {}".format(self.X.shape[0], self.X.shape[1]))
        clf = self.svc.fit(self.X, self.y)
        return clf
    
    #def predict(self, new_y):
    #    self.clf.
        
    def dishonest_cross_validate(self):
        #print("cv over {} x {}".format(self.X.shape[0], self.X.shape[1]))
        scores = cross_val_score(self.svc, self.X, self.y, cv=self.cv_buckets)
        return scores
        
    #def plot(self):
    #    # ONLY 2D !!!
    #    
    #    # create a mesh to plot in
    #    h = .02  # step size in the mesh
    #    x_min, x_max = self.X[:, 0].min() - 1, self.X[:, 0].max() + 1
    #    y_min, y_max = self.X[:, 1].min() - 1, self.X[:, 1].max() + 1
    #    print(x_min, x_max)
    #    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
    #                         np.arange(y_min, y_max, h))
#
    #    # Plot the decision boundary. For that, we will assign a color to each
    #    # point in the mesh [x_min, x_max]x[y_min, y_max].
    #    Z = self.clf.predict(np.c_[xx.ravel(), yy.ravel()])
#
    #    # Put the result into a color plot
    #    Z = Z.reshape(xx.shape)
    #    
    #    fig, ax = plt.subplots(1, 1)
    #    ax.contourf(xx, yy, Z, cmap=plt.cm.coolwarm, alpha=0.8)
#
    #    # Plot also the training points
    #    ax.scatter(self.X[:, 0], self.X[:, 1], c=self.y, cmap=plt.cm.coolwarm)
    #    #ax.set_xlabel('Sepal length')
    #    #ax.set_ylabel('Sepal width')
    #    ax.set_xlim(xx.min(), xx.max())
    #    ax.set_ylim(yy.min(), yy.max())
    #    ax.set_xticks(())
    #    ax.set_yticks(())
    #    ax.set_title(self.title)
#
    #    plt.show()
    #    plt.savefig('../out/{}_classified_{}.png'.format(
    #        self.hypo.batch.tag, self.hypo.batch.onlytag))
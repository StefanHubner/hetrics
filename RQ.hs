module RQ where

import Numeric.LinearAlgebra
import Numeric.GSL.Minimization

type Mat = Matrix Double

data MaxControl = MaxControl { maxit :: Int, epsilon :: Double, huge :: Double, beta :: Double } deriving (Show)
data Model = Model {perc ::Double, endog :: Mat, exog :: Mat} deriving (Show)

defaultMaxCtrl :: MaxControl
defaultMaxCtrl = MaxControl 10 10e-2 10e5 0.97

n :: Int 
n = 1000

mb :: Mat 
mb = trans . fromLists $ [[1, 2, 3]] 

fitfromlist :: [[Double]] -> [[Double]] -> Mat
fitfromlist = curry $ uncurry fitMedianAt0 . uncurry createData

fitMedianAt0 :: Mat -> Mat -> Mat
fitMedianAt0 x y = fit (Model 0.5 y x) (asRow . constant 0 $ rows mb)

createData :: [[Double]] -> [[Double]] -> (Mat, Mat)
createData x u = let (mx, mu) = (g x, g u) in (mx, mx <> mb + mu) 
	where g = trans . (fromLists :: [[Double]] -> Mat)

fromScalar :: Element e => e -> Matrix e 
fromScalar x = asRow . constant x $ 1

lossRQ :: Double -> Mat -> Mat
lossRQ tau u = u * ( t - oneiflt0 u ) 
	where 
		oneiflt0 = mapMatrix (\x -> if x < 0 then 1 else 0)
		t = fromScalar tau

sumLoss :: Double -> Mat -> Double		
sumLoss tau = sumElements . takeColumns 1 . lossRQ tau 

fit :: Model -> Mat -> Mat
fit m x0 = let (_,_,p) = minloss m x0 s0 w0 p0 in p
	where
		w0 = zeros . rows . endog $ m 
		s0 = sumLoss (perc m) $ endog m
		p0 = zeros . (+1) . cols . exog $ m 
		zeros = asColumn . constant 0 

minloss :: Model -> Mat -> Double -> Mat -> Mat -> (Double, Mat, Mat)
minloss m x0 s w par = 	(if s - s' > epsilon defaultMaxCtrl then minloss m x0 else (,,) ) s' w'' par'
	where
		grad = gradient (exog m) x0 par
		(_, w', step') = meketon (Model (perc m) (resid (endog m) grad par) grad) (huge defaultMaxCtrl) w 
		par' = par + (stepsize (Model (perc m) (endog m) grad) par step' `scale` step')
		s' = sumLoss (perc m) . resid (endog m) grad $ par'
		what = predictlin grad $ linearSolveLS grad w'
		(w1,w0) = (maxElement what, -minElement what)
		w'' | w1 > perc m = (perc m / w1) `scale` what  
			| w0 > 1 - perc m = ((1 - perc m) / w0) `scale` what
			| otherwise = what

meketon :: Model -> Double -> Mat -> (Double, Mat, Mat)
meketon m yw w =	if yw - yw' < epsilon defaultMaxCtrl then 
						(yw', w', wbeta)
					else
						meketon m yw' w'
	where	
		d = mapMatrix (\x -> if perc m - x < 0.5 then perc m - x else 1 - perc m + x) w
		wbeta = linearSolveLS (exog m * d) (endog m * d)
		wresid = resid (endog m) (exog m) wbeta
		s =  wresid * d * d 
		alpha = let t = fromScalar . perc $ m in maxElement $ liftMatrix2 max (s/(t - w)) (-s/(1 - t + w))
		w' = w + ((beta defaultMaxCtrl/) . max (epsilon defaultMaxCtrl) $ alpha) `scale` s 
		yw' = sumLoss (perc m) wresid 

resid :: Mat -> Mat -> Mat -> Mat
resid y x par = y - predictlin x par 

predictlin :: Mat -> Mat -> Mat
predictlin x par = x <> par 

gradient :: Mat -> Mat -> Mat -> Mat
gradient x x0 _ = fromBlocks [[asColumn . constant 1 . rows $ x, x - x0]] 

stepsize :: Model -> Mat -> Mat -> Double					
stepsize m par step' = let (sol, _) = minimize NMSimplex (epsilon defaultMaxCtrl) (maxit defaultMaxCtrl) [1.0] obj [0.5] in head sol 
	where
		obj :: [Double] -> Double
		obj l = sumLoss (perc m) $ endog m - predictlin (exog m) (par + head l `scale` step') 
-- uniMinimize BrentMini (epsilon defaultMaxCtrl) (maxit defaultMaxCtrl) obj 0.5 0.0 10.0 in sol 

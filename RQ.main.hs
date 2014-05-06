import Control.Monad
import Numeric.LinearAlgebra
import System.Random.MWC
import System.Random.MWC.Distributions
import RQ 

main :: IO ()
main = do
	s <- create
	u:x <- replicateM (rows mb + 1) . replicateM n . normal 0 1 $ s 
	print $ toList $ flatten $ fitfromlist x [u] 
	return ()

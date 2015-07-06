--Position of the point of view
DECLARE @x FLOAT
DECLARE @y FLOAT
DECLARE @z FLOAT
--Snapshot to be analized
DECLARE @snapnum INT
--Hubble constant
DECLARE @H FLOAT
-- Little h
DECLARE @lh float
--Neighbour
DECLARE @neigh INT
--Limits the box of halos
DECLARE @limx FLOAT
DECLARE @limy FLOAT
DECLARE @limz FLOAT
DECLARE @limMx FLOAT
DECLARE @limMy FLOAT
DECLARE @limMz FLOAT


--Asignment of variables
SET @lh = 0.7
SET @snapnum = 63
SET @neigh = 4
SET @x = 0
SET @y = 0
SET @z = 0
SET @limx = 16
SET @limy = 16
SET @limz = 16
SET @limMx = 40
SET @limMy = 40
SET @limMz = 40
SET @H = 100

DECLARE @hid INT
DECLARE @n INT
SET @hid = 1;
SET @n = 1;



SELECT *
FROM
      (SELECT *, ROW_NUMBER() OVER ( PARTITION BY f.haloID ORDER BY f.SIGMA ASC ) AS RN
      FROM
             (SELECT (
                        POWER(
                              (
                                 (D.x - @x)*(M.x - @x)
                               + (D.y - @y)*(M.y - @y)
                               + (D.z - @z)*(M.z - @z)

                             )*(M.x - @x)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                               - (D.x - @x)
                               , 2 )
                               +
                        POWER(
                              (
                                 (D.x - @x)*(M.x - @x)
                               + (D.y - @y)*(M.y - @y)
                               + (D.z - @z)*(M.z - @z)

                             )*(M.y - @y)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                               - (D.y - @y)
                               , 2 )
                               +
                        POWER(
                              (
                                 (D.x - @x)*(M.x - @x)
                               + (D.y - @y)*(M.y - @y)
                               + (D.z - @z)*(M.z - @z)

                             )*(M.z - @z)/( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )
                               - (D.z - @z)
                               , 2 )

                        ) AS SIGMA,
                        M.np AS NP,
                        M.haloID AS haloID
              FROM DeLucia2006a D,
                   MPAHalo M
              WHERE --Selects snapshot
                    D.snapnum = @snapnum
                    AND M.snapNum = @snapnum
                    -- Limits the box of haloes
                    AND M.x > @limx AND  M.y > @limy AND M.z > @limz
                    AND M.x < @limMx AND  M.y < @limMy AND M.z < @limMz
                    --Omits galaxies with r abs magnitude < -19
                    AND D.mag_r < -19
                    --Selects only galaxies that fulfill the velocity cut
                    AND ABS(
                         ((@H*(D.x - M.x)/@lh + D.velX - M.velX)*(M.x - @x)+
                         (@H*(D.y - M.y)/@lh + D.velY - M.velY)*(M.y - @y) +
                         (@H*(D.z - M.z)/@lh + D.velZ - M.velZ)*(M.z - @z))/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) )

                           )  < 500
              ) f
       )g
WHERE RN = @neigh

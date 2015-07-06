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
      (SELECT  f.DIST,
               4.0/(DIST*PI()) AS SIGMA,
               f.NP,
               ROW_NUMBER() OVER ( PARTITION BY f.haloID ORDER BY f.DIST ASC ) AS RN
      FROM
            (SELECT
                  POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dx - (c.gx - @x)  , 2 )  +
                  POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dy - (c.gy - @y)  , 2 )  +
                  POWER(  ((c.gx - @x)*c.dx + (c.gy - @y)*c.dy + (c.gz - @z)*c.dz)*c.dz - (c.gz - @z)  , 2 )
                  AS DIST,
                  c.haloID AS haloID,
                  c.NP AS NP
            FROM
                -- Selects the position, velocities and position, velocities of galaxies and haloes and the unit vectors to haloes
                (SELECT M.haloID AS haloID,
                        M.np AS NP,
                        D.x/@lh AS gx,
                        D.y/@lh AS gy,
                        D.z/@lh AS gz,
                        M.x/@lh AS hx,
                        M.y/@lh AS hy,
                        M.z/@lh AS hz,
                        --Difference in peculiar velocities
                        D.velX - M.velX AS dvx,
                        D.velY - M.velY AS dvy,
                        D.velZ - M.velZ AS dvz,
                        (M.x - @x)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dx,
                        (M.y - @y)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dy,
                        (M.z - @z)/SQRT( POWER(M.x - @x,2) + POWER(M.y - @y,2) + POWER(M.z - @z,2) ) AS dz
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
                ) c
            WHERE --Only counts galxies with relative radial velocities to the halo < 500 Km/s
                 ABS(
                     -- Radial component of the difference in peculiar velocities
                     c.dvx*c.dx + c.dvy*c.dy + c.dvz*c.dz +
                     -- Radial component of the diference in velocities due to Hubbles law
                     @H*((c.gx - c.hx)*dx + (c.gy - c.hy)*dy + (c.gz - c.hz)*dz)

                   )  < 500
           ) f
      ) g
WHERE RN = @neigh

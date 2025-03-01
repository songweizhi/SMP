# https://www.rdocumentation.org/packages/vegan/versions/1.11-0/topics/betadisper

library('vegan')

data(varespec)
varespec

## Bray-Curtis distances between samples
dis <- vegdist(varespec)
dis

## First 16 sites grazed, remaining 8 sites ungrazed
groups <- factor(c(rep(1,16), rep(2,8)), labels = c("grazed","ungrazed"))
groups

## Calculate multivariate dispersions
mod <- betadisper(dis, groups)
mod

## Perform test
anova(mod)

## Plot the groups and distances to centroids on the first two PCoA axes
plot(mod)

## Draw a boxplot of the distances to centroid for each group
par(cex.axis=1, mar=c(5, 5, 5, 5))
boxplot(mod,
        # col = c("lightgreen", "orange"),
        xlab = 'Distance to group centroids',
        ylab = '', # Sponge species (environments) 
        las=1,  # control direction of group label
        horizontal=TRUE)

#' title: "Crop Analysis Q3 2013"
#' author: "John Smith"
#' date: "May 3rd, 2014"
#' head1
#'   head2
#'      head3
#'         head4
#' 
#' 
#' 
# Section 1
# Class 5
read.table("bggn213_05_rstats/weight_chart.txt")
baby <- read.table("bggn213_05_rstats/weight_chart.txt", header = TRUE)
#View(baby)
plot(baby)
plot(baby, type="b", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Some title")

#section 1B
feat <- read.table("bggn213_05_rstats/feature_counts.txt", sep="\t", header = TRUE)
#View( feat )
barplot( feat$Count, names.arg=feat$Feature, main = "some title", horiz = TRUE, las=1  )

# Section C
# hist(10000, 0)
rnorm(10000, 0)
# hist(rnorm(10000, 0))
read.delim("bggn213_05_rstats/male_female_counts.txt")
mfcount <- read.delim("bggn213_05_rstats/male_female_counts.txt") 
barplot(mfcount$Count, col=rainbow(10))

nrow(mfcount)

# 3
expr <- read.delim("bggn213_05_rstats/up_down_expression.txt", header = TRUE)
plot(expr$Condition1, expr$Condition2, col=expr$State)
# how may genes went up and down arround?
table(expr$State)
plot(1:5, pch=15, col=2)
?palete
??palette

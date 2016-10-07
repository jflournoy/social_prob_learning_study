install.packages('ndl')
library(ndl)

data(lexample)

library(tidyr)
library(dplyr)

somedata <- read.csv('~/code/socprobmisc/split-130-1468342589910.csv')

justthefacts <- somedata %>% filter(stimulus != '') %>%
	extract(stimulus, c('image_file', 'left', 'right'), '.*([fm]_f42887.*png).*left>(\\w+)</div>.*right>(\\w+)</div>') %>% 
	select(image_file, left, right, optimal_response, correct_response)

table(justthefacts$image_file, justthefacts$correct_response)


#Well, that's exciting. SON OF A
data.frame(table(justthefacts$image_file, justthefacts$correct_response)) %>%
	spread(Var2, Freq) %>%
	rowwise() %>%
	mutate(prop37 = `37`/(`37`+`39`),
	       prop39 = `39`/(`37`+`39`))

orthoCoding(lexample$Word, grams=1)

splitCues <- data.frame(Cues=rep(c('h_t_z', 'h_t_y', 
				   'd_l_x', 'd_l_w', 
				   'p_n_v', 'p_n_u'), each=2),
			Outcomes=c(rep(c('hungry','thirsty'), 2),
				   rep(c('dating','looking'), 2),
				   rep(c('popular', 'unpopular'), 2)),
			Frequency=rep(c(51, 13, 13, 51), 3),
			stringsAsFactors = F)

rw.test <- RescorlaWagner(splitCues, nruns=1, traceCue='z', traceOutcome='thirsty')
plot(rw.test)

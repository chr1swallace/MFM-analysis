#!/usr/bin/ruby
require 'optparse'
require 'fileutils'
require '/home/cew54/slurmer/Qsub.rb'
require 'pp'
require 'date'

load("~/DIRS.txt")
#ROOT="/home/cew54/scratch"
#IMPUTEDIR="/home/cew54/scratch/FM-impute"
QCIMPUTEDIR="#{ROOT}/FM-impute_qc"
#GUESSDIR="#{ROOT}/FM-GUESS"
#SINGLESNPDIR="/home/cew54/scratch/FM-ss"
#OUTPUTDIR="/home/cew54/scratch/FM/output-new"

OPTIONS = {}
OPTIONS[:int] = false
OPTIONS[:run] = false
OPTIONS[:all] = false
OptionParser.new do |opts|
  opts.banner = "Usage: runguess.rb [OPTIONS] COMMAND [DIR]"

  opts.on("-a", "--[no-]all", "(Re-)run all outputs regardless of whether it already exists") do |i|
    OPTIONS[:all] = i
  end
  
  opts.on("-i", "--[no-]interactive", "Run each job in serial on the interactive log in node") do |i|
    OPTIONS[:int] = i
  end
  
  opts.on("-r", "--autoRun", "Run each job on the queue without user input") do |i|
    OPTIONS[:autorun] = i
  end
  
  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    OPTIONS[:verbose] = v
  end

  opts.on("-h", "--help", "Show help") do |h|
    OPTIONS[:help] = h
  end
  
  # opts.on("-q", "--[no-]queue", "Use qCom.sh in output") do |v|
  #   OPTIONS[:queue] = v
  # end
  # opts.on("-n", "--nohup", "Prepend with nohup") do |v|
  #   OPTIONS[:nohup] = v
  # end
  # opts.on("-b", "--[no-]background", "Append with &") do |v|
  #   OPTIONS[:background] = v
  # end
end.parse!
COMMAND = ARGV.shift

def usage()
  puts OPTIONS
  puts "Usage: runguess.rb [OPTIONS] COMMAND [DIR]

  COMMANDS are:
      (run things)
           prep :: prepare GUESS input
           guess :: run guess
           expand :: expand guess output
           qc :: generate qc summaries
           group :: generate group.RData
           plot :: generate summx and summary results plots
           MTFM :: run MTFM
           plot2 :: generate summx and summary results plots including for MTFM
           
           mann :: single-snp-pvalues.csv
           conditional :: conditional.RData
           condfm :: conditional-finemap.RData
           nsnp :: nsnp.pdf
           
           collate :: gather output from plot2 in ~/scratch/cvs-table
           html :: format for html output ~/scratch/MFM-output

      (manage things)
           missing :: report directories with guess results but no snpmod (expand fail)
           guesscountbad :: count incomplete guess runs
           guessrmbad :: remove incomplete guess runs

           clean :: remove slurm files      
           summary :: describe what is done/undone
      
      (old things for CVS paper)
           tables :: make parts of supplementary tables
           prepint :: prepare international sample data for same snps fed to GUESS
           expandall :: expand for *all* diseases not just 4 with most data (output is allmod.RData)
           fixprior :: fix priors and recompute SM2
           fixpriorall :: fixprior on these 
           group :: generate alt-group.RData
           plotgroup :: compare stepwise & GUESSFM groups in plot format 
           roundup :: summarise full results to date
           haplotypes :: generate haplotype analyses

  "
end
if OPTIONS[:help] then
  usage()
  exit 0
end

comfile="runguess-#{COMMAND}.sh"
args = {:job=>COMMAND,
        :tasks=>'1',
        :cpus=>'1',
           :autorun=>OPTIONS[:autorun],
           :excl=>" "}

# q=Qsub.new("runguess-#{command}.sh",
#            :job=>command,
#            :tasks=>'16',
#            :time=>options["t"],
#            :account=>options["a"],
#            :excl=>options["x"],
#            :autorun=>options["r"])


def command_prep(region)  
  logf = "log/prepguess-#{region}.log"
  inf = "#{QCIMPUTEDIR}/imputed-#{region}.RData"
  "./R/prepguess.R --args file=#{inf} > #{logf} 2>&1"
end
def command_prepint(region)  
  logf = "log/prepint-#{region}.log"
  inf = "#{QCIMPUTEDIR}/imputed-#{region}.RData"
  "./R/prep-intdata.R --args file=#{inf} > #{logf} 2>&1"
end

def command_qc(f)  
  "./R/GUESS-qc.R --args d=#{f}"
end

def command_guess(f)
  "/home/cew54/R/x86_64-pc-linux-gnu-library/3.3/R2GUESS/bin/GUESS -history -X #{f}/X_50000 -Y #{f}/Y_50000 -nsweep 80000 -burn_in 4545 -out #{f}/out_50000 -par #{f}/par.xml -top 5000 -init #{f}/init_50000 -Egam 2 -Sgam 1.41013212943979 -n_chain 5 > #{f}/log 2>&1"
end
def command_expand(f)
  region=get_region(f)
  "R/GUESS-expand-models.R --args d=#{f} > log/GUESS-expand-models-#{region}.Rout 2>&1"
end
def command_fixprior(f)
  region=get_region(f)
  "R/GUESSFM-fixpriors.R --args d=#{f} > log/GUESS-fixpriors-#{region}.Rout 2>&1"
end
def command_plot(d)
  region=get_region(d)
  "R/GUESS-resultplots.R --args d=#{d} > log/GUESS-resultplots-#{region}.Rout 2>&1"
end
def command_haplotype(d)
  region=get_region(d)
  "R/GUESS-haplotypes.R --args d=#{d} > log/GUESS-haplotypes-#{region}.Rout 2>&1"
end
def command_conditional(d)
  region=get_region(d)
  "R/conditional.R --args d=#{d} > log/conditional-#{region}.Rout 2>&1"
end
def command_R(region,com,args='',parent='')
  if(parent=='..') then
    "R/#{com}.R --args d=#{CVSROOT}/#{region} #{args} > log/#{com}-#{region}.Rout 2>&1"
  else 
    "R/#{com}.R --args d=#{CVSROOT}/#{region}/GUESS #{args} > log/#{com}-#{region}.Rout 2>&1"
  end
end


################################################################################

## existing files
def cpifnewer(f,d)
  of = d + "/" + File.basename(f)
  if(File.exists?(f)) then
    FileUtils.cp f, d unless  File.exists?(of)
    FileUtils.cp f, d if FileUtils.uptodate?(f,[of]) #or FileUtils.cp f, d
  end
end

def listfiles(dir,patt)
  files = Dir.glob(dir + '/' + patt)
  # puts "FILES in #{dir}" if OPTIONS[:verbose]
  # puts files if OPTIONS[:verbose]
  return files
end
def listbadguess()
  guessdirs = listfiles(GUESSDIR, '[0-9]*')
  delfiles= []
  guessdirs.each { |d|
    guessout = listfiles(d, '[iAGCIRTJM]*/out_50000_output_best_visited_models.txt')
    guessout = guessout.select{ |f| `wc -l "#{f}"`.strip.split(' ')[0].to_i == 1  }   
    guessout_dirs = guessout.map{ |f| File.dirname(f) }
    delfiles.push(*guessout_dirs) if guessout_dirs.length() > 0
  }
  delfiles
end
def getdirs(files)
  files.map { |f|
    File.dirname(f)
  }
end
def shortensummx(d) 
    m = d.match /(.+\d+[pq]-\d+-\d+)-/
    m ? m[1] : d
end
def get_region(s)
  /(\d+[pq]-\d+-[\d\+rs]+)/.match(s)[1]
end
################################################################################
guessdirs = listfiles(CVSROOT, '[0-9]*/GUESS').select { |d|
  guessout = listfiles(d, '[iAGCIRTJM]*/out_50000_features.txt')
  guessout = guessout.select{ |f| `wc -l "#{f}"`.strip.split(' ')[0].to_i > 1  }
  guessout = guessout.sort_by{ |f| File.mtime(f) }.reverse!
  skipfile = d + '/skip'
  guessout.length()>0 && !(File.exists?(skipfile) && (File.mtime(skipfile) > File.mtime(guessout[0])))
}

## command
commands = []

case COMMAND
when "clean"
  rmfiles=listfiles(".", "runguess-*.sh*") +
          listfiles(".", "slurm-*.sh*") +
          listfiles(".", "machine.file.*")
  if(rmfiles.length() > 0) then
    puts "removing #{rmfiles.length()} files"
    puts(rmfiles)
    rmfiles.map { |f| File::delete(f) }
  end

when "prep"
  puts "--- prepguess.R ---"
  args[:time] = "1:00:00"
  args[:cpus] = '12' # memory issues?
  impfiles = listfiles(IMPROOT, 'imputed-*.RData')
  guessin = listfiles(CVSROOT, '*/GUESS').map {|f| File.dirname(f)}
  todo = impfiles.map { |f| get_region(f) } - guessin.map { |f| get_region(f) } # no GUESS DIR created
  # guessout = listfiles(CVSROOT, '*/GUESS/[GCIRTJM]*/out_50000_features.txt').map {|f| File.dirname(f)}
  # guessskip = listfiles(CVSROOT, '*/GUESS/[GCIRTJM]*/out_50000_features.txt').map {|f| File.dirname(f)}
  todo2 = guessin.select{ |d|
    done=listfiles(d,'GUESS/[ACRTJM]*/out_50000_features.txt').length
    skipped=listfiles(d,'GUESS/[ACRTJM]*-skip').length
    done + skipped < 6
  }.map { |f| get_region(f) }
  # commands = todo.map { |r| command_prep(r) }
  commands = (todo + todo2).map { |r| command_prep(r) }

when "prepint"
  puts "--- prepint.R ---"
  args[:cpus] = '1'
  guessin = listfiles(CVSROOT, '*/GUESS/all-data.RData').map {|f| File.dirname(File.dirname(f)) }
  intin = listfiles(CVSROOT, '*/GUESS/int-data.RData').map {|f| File.dirname(File.dirname(f)) }
  todo = guessin.map { |f| get_region(f) } - intin.map { |f| get_region(f) }
  commands = todo.map { |r| command_prepint(r) }

when "mann"
  puts "--- manhattan data ---"
#  args[:excl] = " --exclusive " # NB - ICOELIAC and RA fail on shared nodes due to memory
  guessin = listfiles(CVSROOT, '*/GUESS/all-data.RData').map {|f| File.dirname(File.dirname(f)) }
  guessout = listfiles(CVSROOT, '*/single-snp-pvalues.csv').map {|f| File.dirname(f)}
  guesstodo = guessin - guessout
  commands = guesstodo.map { |d| command_R(get_region(d),"manhattans") }

when "guess"
  puts "--- GUESS ---"
#  args[:excl] = " --exclusive " # NB - ICOELIAC and RA fail on shared nodes due to memory
  args[:time] = "48:00:00"
  args[:tasks] = '1'
  args[:cpus] = '4' # bigger memory
  guessin = listfiles(CVSROOT, '*/GUESS/[iAGCIRTJM]*/init_50000').map {|f| File.dirname(f)}
  guessout = listfiles(CVSROOT, '*/GUESS/[iAGCIRTJM]*/out_50000_features.txt').map {|f| File.dirname(f)}
  guesstodo = guessin - guessout
  if ARGV.length() > 0 then
    guesstodo = guesstodo & ARGV
  end
  commands = guesstodo.map {|f| command_guess(f) }

when "qc"
  puts " --- GUESS-qc.R ---"
   args[:cpus] = 1
   args[:tasks] = 1
  args[:time] = "00:50:00"
  commands = guessdirs.map { |d| command_R(get_region(d),"GUESS-qc") }

when "expand"
  # args[:excl] = " --exclusive " # NB - now paralellise within R
  args[:time] = "24:00:00"
  args[:cpus] = "12"
  # guessout = listfiles(GUESSDIR, '*/[GCIRTJM]*/out_50000_features.txt').map {|f| File.dirname(f)}
  # expandout = listfiles(GUESSDIR, '*/snpmod.RData').map {|f| File.dirname(f)}
## run when no expandout for a guess input dir, or when the expandout is older than any corresponding guessout and the guessout is finished
  guessdirs.each { |d|
    snpmod = listfiles(d,"../snpmod-99.RData")
    guessout = listfiles(d, '[iAGCIRTJM]*/out_50000_features.txt')
    guessout = guessout.select{ |f| `wc -l "#{f}"`.strip.split(' ')[0].to_i > 1  }
    guessout = guessout.sort_by{ |f| File.mtime(f) }.reverse!
    if snpmod.length()==0 || File.mtime(snpmod[0]) < File.mtime(guessout[0]) then
      commands.push( command_R(get_region(d), "expand") )
      puts get_region(d) if OPTIONS[:verbose]
    end
  }

when "group"
  args[:time] = "0:45:00"
  args[:tasks] = '1'
  args[:cpus] = 12
  summx = listfiles(CVSROOT, '[0-9]*/snpmod-99.RData') #.map { |f| get_region(f) }
  output = listfiles(CVSROOT, "[0-9]*/snpmod-99-groups.RData").map {|f| f.sub("-groups","") } #get_region(f) }
  if(OPTIONS[:all])
    todo = summx
  else
    todo = summx - output
  end
  todor = todo.map { |f| get_region(f) }
  thr = todo.map { |f| /snpmod-([0-9]+).RData/.match(f)[1] }
  puts "group todo " + todo.length.to_s
  todor.each_with_index {|r,i|
    commands.push(command_R(r,"group-snps-v3","cpp.thr=#{thr[i]}",parent=".."))
  }

when "MTFM"
  args[:time] = "6:00:00"
  args[:tasks] = '1'
  args[:cpus] = 3 # parallise over 3 values for kappa
  input = listfiles(CVSROOT, '[0-9]*/snpmod-99-groups.RData').map { |f| get_region(f) }
  input2 = listfiles(CVSROOT, '[0-9]*/conditional.RData').map { |f| get_region(f) }
  if (input - input2).length() > 0
    puts "!!! run ./run.rb conditional"
    puts (input - input2)
  end
  input=input & input2
  skips = listfiles(CVSROOT, '[0-9]*/skip-mtfm').map { |f| get_region(f) }
  output = listfiles(CVSROOT, "[0-9]*/MTFM.RData").map {|f| get_region(f) }
  output2 = listfiles(CVSROOT, "[0-9]*/iMTFM.RData").map {|f| get_region(f) }
  todo = input - skips
  if(!OPTIONS[:all])
    todo = todo - output - output2
  end
  puts "MTFM todo " + todo.length.to_s
  todo.each {|r|
    commands.push(command_R(r,"MTFM",'',parent=".."))
  }

when "nsnp"
  commands.push("Rscript R/nsnps.R")
  
when "plot"
  puts " --- R/MTFM-plots.R --- "
  ## run locally
  args[:tasks] = '1'
  # args[:time] = '04:00:00'
  snpmod = listfiles(CVSROOT, '*/snpmod-99.RData')
  summx = listfiles(CVSROOT, '*/plots/summx2-best.pdf')
  snpmodr = snpmod.map {|f|
    get_region(f)
  }
  summxr = summx.map { |f|
    get_region(f)
  }
  if(OPTIONS[:all])
    todo = snpmodr
  else 
    ## easy - snpmods with no matching summx
    nomatch = snpmodr - summxr
    ## for matches, check if summx older than snpmod
    matches = snpmodr & summxr
    matches = matches.select { |r|
    snpmod_t = File.mtime(CVSROOT + "/" + r + "/snpmod-99.RData")
    summx_f = summx.find { |f| /#{Regexp.quote(r)}/ =~ f }
    summx_t = File.mtime(summx_f)
    snpmod_t > summx_t
  }
  todo = nomatch + matches
  end
  puts "running " + todo.length.to_s + " regions"
  todo.each {|r|
    comm=command_R(r, "MTFM-plots")
    commands.push(comm)
    #puts comm
    #system("#{comm}")
  }
  
when "plot2"
  puts " --- R/MTFM-plot-part2.R --- "
  ## run locally
  args[:tasks] = '1'
  args[:cpu] = '6'
  snpmod = listfiles(CVSROOT, '*/snpmod-99.RData')
  summx = listfiles(CVSROOT, '*/plots/haplotypes.pdf')
  mtfm = listfiles(CVSROOT, '*/MTFM.RData')
  snpmodr = snpmod.map {|f|
    get_region(f)
  }
  summxr = summx.map { |f|
    get_region(f)
  }
  if(OPTIONS[:all])
    todo = snpmodr
  else 
    ## easy - snpmods with no matching summx
    nomatch = snpmodr - summxr
    ## for matches, check if summx older than snpmod
    matches = snpmodr & summxr
    matches = matches.select { |r|
    snpmod_t = File.mtime(CVSROOT + "/" + r + "/snpmod-99.RData")
    summx_f = summx.find { |f| /#{Regexp.quote(r)}/ =~ f }
    summx_t = File.mtime(summx_f)
    snpmod_t > summx_t
  }
  todo = nomatch + matches
  end
  puts "running " + todo.length.to_s + " regions"
  puts "after this, run"
  puts "grep variability log/MTFM-plot-part2-*.Rout"
  todo.each {|r|
    comm=command_R(r, "MTFM-plot-part2")
    commands.push(comm)
    #puts comm
    #system("#{comm}")
  }
  
when "tables"
  args[:tasks] = '1'
  summx = listfiles(CVSROOT, '*/plots/summx2-best.pdf')
  summxr = summx.map { |f|
    get_region(f)
  }
  tabs = listfiles("/home/cew54/scratch/cvs-table","[0-9]*[pq]-*.tex")
  tabr = tabs.map { |f|
    get_region(f)
  }
  todo=summxr
  if(!OPTIONS[:all])
    todo = summxr - tabr
  end
  puts "running supplementary-table.R for " + todo.length.to_s + " regions"
  todo.each {|r|
    comm=command_R(r, "supplementary-table")
    commands.push(comm)
    puts(r) if OPTIONS[:verbose]
    #puts comm
    #system("#{comm}")
  }
  
when "html"
  todo = listfiles(CVSROOT, '*/plots/haplotypes.pdf').map { |f|
    get_region(f)
  }
    todo = todo - [ '10p-6030243-6169685' ]
todo = todo.sort do |a, b|
    af=a.sub(/p|q/,"").sub(/-/,".").to_f # only first - gets switched to ., then to_f will drop everything after the remaining -
    bf=b.sub(/p|q/,"").sub(/-/,".").to_f
    af <=> bf
  end
  todo = todo - %w(2q-204446380-204816382+rs117701653 10p-6030243-6169685)
 
  todo.each do |r|
      puts r
      od ="/home/cew54/scratch/MFM-output/#{r}"
      Dir.mkdir(od) unless File.exists?(od)
      ## cp files over if needed
      # %w(tables.tex groupld.pdf summx-mpp.pdf summx-impp.pdf summx-ind.pdf haplotypes.pdf).each { |f|
      %w(groupld.pdf summx-mpp.pdf summx-impp.pdf summx-ind.pdf haplotypes.pdf nsnps.pdf).each { |f|
        cpifnewer "#{CVSROOT}/#{r}/plots/#{f}",od
      }
      # cpifnewer "#{CVSROOT}/#{r}/GUESS/qc3-plots.pdf","#{od}/nsnps.pdf"
  end

  # system("cp /home/cew54/scratch/MFM-output/2q-204446380-204816382+rs117701653/* /home/cew54/scratch/MFM-output/2q-204446380-204816382")
  system("Rscript ./R/html.R")
  puts "now "
  puts "cd /home/cew54/scratch/MFM-output"
  puts "git commit -a -m \"update #{Date.today}\""
  puts "git push"
#   end 
#   FT.close()
#   FF.close()
#   FH.close()
#   FA.close()

  

when "collate"
  todo = listfiles(CVSROOT, '*/plots/haplotypes.pdf').map { |f|
    get_region(f)
  }
  todo = todo - [ '10p-6030243-6169685' ]
  todo = todo.sort do |a, b|
    af=a.sub(/p|q/,"").sub(/-/,".").to_f # only first - gets switched to ., then to_f will drop everything after the remaining -
    bf=b.sub(/p|q/,"").sub(/-/,".").to_f
    af <=> bf
  end
  ## create list of good regions
  File.open("/home/cew54/scratch/cvs-table/regions.txt", 'w') {|f|
    todo.each do |r|
      f.puts r
    end
  }

  FT=File.open("/home/cew54/scratch/cvs-table/index-tables.tex", 'w')
  FH=File.open("/home/cew54/scratch/cvs-table/index-haps.tex", 'w')
  FA=File.open("/home/cew54/scratch/cvs-table/index.tex", 'w')
  FF=File.open("/home/cew54/scratch/cvs-table/index-figures.tex", 'w')
  todo.each do |r|
      puts r
      od ="/home/cew54/scratch/cvs-table/#{r}"
      Dir.mkdir(od) unless File.exists?(od)
      ## cp files over if needed
      %w(tables.tex groupld.pdf summx-mpp.pdf summx-impp.pdf summx-ind.pdf haplotypes.pdf).each { |f|
        cpifnewer "#{CVSROOT}/#{r}/plots/#{f}",od
      }

      ## add lines to .tex file
      ## tables
      FT.puts("\\clearpage
\\markright{#{r}}
\\centerline{\\textbf{#{r}}}
\\input{cvs-table/#{r}/tables.tex}
")
      FA.puts("
\\clearpage
\\markright{#{r}}
\\centerline{\\textbf{#{r}}}
\\input{cvs-table/#{r}/tables.tex}
")
      FF.puts("
\\clearpage
\\markright{#{r}}
")
      FH.puts("
\\clearpage
\\markright{#{r}}
")
       ## LD
      if File.exists?("#{CVSROOT}/#{r}/plots/groupld.pdf") then
        FT.puts("\\centerline{\\includegraphics{cvs-table/#{r}/groupld.pdf}}
\\clearpage
")
      FA.puts("\\centerline{\\includegraphics{cvs-table/#{r}/groupld.pdf}}
\\clearpage
") 
      end

      ## haps + first summx
      FA.puts("\\clearpage
Haplotypes\\newline
\\includegraphics[height=0.9\\textheight]{cvs-table/#{r}/haplotypes.pdf}\\newline
\\clearpage
Independent GUESSFM analysis\\newline
\\includegraphics[angle=90,height=0.95\\textheight]{cvs-table/#{r}/summx-ind.pdf}\\newline
")
      FH.puts("\\clearpage
\\includegraphics[height=0.9\\textheight]{cvs-table/#{r}/haplotypes.pdf}\\newline
")
FF.puts("\\clearpage
Independent GUESSFM analysis\\newline
\\includegraphics[angle=90,height=0.95\\textheight]{cvs-table/#{r}/summx-ind.pdf}\\newline
")
      
       ## optional MTFM bit
      if File.exists?("#{CVSROOT}/#{r}/plots/summx-mpp.pdf") then
           FA.puts("\\clearpage
Joint MTFM UK analysis\\newline
\\includegraphics[angle=90,height=0.95\\textheight]{cvs-table/#{r}/summx-mpp.pdf}\\newline
") 
           FF.puts("\\clearpage
Joint MTFM UK analysis\\newline
\\includegraphics[angle=90,height=0.95\\textheight]{cvs-table/#{r}/summx-mpp.pdf}\\newline
") 
      end

      ## optional international bit
      if File.exists?("#{CVSROOT}/#{r}/plots/summx-impp.pdf") then
          FA.puts("\\clearpage
Joint MTFM International analysis\\newline
\\includegraphics[angle=90,height=0.95\\textheight]{cvs-table/#{r}/summx-impp.pdf}\\newline
") if File.exists?("#{CVSROOT}/#{r}/plots/summx-impp.pdf")
          FF.puts("\\clearpage
Joint MTFM International analysis\\newline
\\includegraphics[angle=90,height=0.95\\textheight]{cvs-table/#{r}/summx-impp.pdf}\\newline
") if File.exists?("#{CVSROOT}/#{r}/plots/summx-impp.pdf")
      end
  end 
  FT.close()
  FF.close()
  FH.close()
  FA.close()

when "condfm"
  args[:tasks] = '1'
  args[:cpus] = '4'
  args[:time] = '08:00:00'
  # snpmod = listfiles(CVSROOT, '{[1-9],[1-9][0-9]}[pq]-*').map {|f| get_region(f) }
  todo = listfiles(CVSROOT, '*/conditional.RData').map {|f| get_region(f) }
  done = listfiles(CVSROOT, '*/conditional-finemap.RData').map {|f| get_region(f) }
  if(!OPTIONS[:all])
    todo = todo - done
  end
  todo.each { |r|
    commands.push(command_R(r, "conditional-finemap"))
    puts r if OPTIONS[:verbose]
  }
  


when "conditional"
  args[:tasks] = '1'
  args[:cpus] = '1'
  args[:time] = '04:00:00'
  # snpmod = listfiles(CVSROOT, '{[1-9],[1-9][0-9]}[pq]-*').map {|f| get_region(f) }
  todo = listfiles(CVSROOT, '*/GUESS').map {|f| get_region(f) } #File.dirname(File.dirname(f))) }
  cond = listfiles(CVSROOT, '*/conditional.RData').map {|f| get_region(f) }
  if(!OPTIONS[:all])
    todo = todo - cond
  end
  todo.each { |r|
    commands.push(command_R(r, "conditional-v2"))
    puts r if OPTIONS[:verbose]
  }
  

# when "tables"
  
when "haplotypes"
  args[:tasks] = '1'
  todo = listfiles(CVSROOT, '*/plots/summw2.csv').map { |f| get_region(f) }
  haplotypes = listfiles(CVSROOT, '*/plots/haplotypes.pdf').map {|f| get_region(f) }
  chaplotypes = listfiles(CVSROOT, '*/cond-haplotypes.pdf').map {|f| get_region(f) }
  if(!OPTIONS[:all])
    todo = todo - haplotypes
  end
  puts "haplotypes todo " + todo.length.to_s
  todo.each {|r|
    commands.push(command_R(r,"GUESS-haplotypes"))
  }
  # todo = summx - chaplotypes
  #   puts "conditional haplotypes todo " + todo.length.to_s
  #   todo.each {|r|
  #     commands.push(command_R(r,"GUESS-haplotypes","cond"))
  #   }

when "summary"
  rdirs = listfiles(CVSROOT, '[0-9]*')
  puts "region count: " + rdirs.length().to_s()
  cond = listfiles(CVSROOT, "*/conditional.RData")
  puts "conditional: " + cond.length().to_s()
  parfiles = listfiles(CVSROOT, '[0-9]*/GUESS/[iAGCIRTJM]*/par.xml')
  pardirs = parfiles.map {|f| File.dirname( File.dirname( File.dirname(f) ) ) }.uniq
  puts "GUESS runs prepared: " + parfiles.length().to_s() + " (#{pardirs.length} regions)"
  skipdirs = listfiles(CVSROOT,  '[0-9]*/GUESS/skip').map { |f|
    File.dirname( File.dirname(f) )
  } & pardirs
  
#  guessdirs = guessdirs - skipped.map { |s| s.sub("/skip","") }
  ngood=0
  nbad=0
  nrgood=0
  nrbad=0
  pardirs.each { |d|
    guessin = listfiles(d, 'GUESS/[iAGCIRTJM]*/par.xml')
    guessout = listfiles(d, 'GUESS/[iAGCIRTJM]*/out_50000_features.txt')
    guessout_good = guessout.select{ |f| `wc -l "#{f}"`.strip.split(' ')[0].to_i > 1  }
    ni = guessin.length()
    ng =  guessout_good.length()
    nb = guessin.length() - guessout_good.length()
    ngood += ng
    nrgood += 1 if ng > 0
    nbad += nb
    nrbad += 1 if nb > 0
  }
  puts "   complete guess runs: #{ngood} (#{nrgood} regions)"
  puts "   incomplete guess runs: #{nbad} (#{nrbad} regions)"
  puts "   regions skipped after GUESS: " + skipdirs.length().to_s()
  puts "   regions in guessdirs to be processed: " + guessdirs.length().to_s()
  snpmod = listfiles(CVSROOT, '*/snpmod-99.RData')
  groups = listfiles(CVSROOT, '*/snpmod-99-groups.RData')
  summx = listfiles(CVSROOT, '*/summx2-99.RData')
  # skipped = listfiles(CVSROOT, '*/GUESS/skip')
  mtfm = listfiles(CVSROOT, "*/MTFM.RData")
  mtfmskip = listfiles(CVSROOT, "*/skip-mtfm")
  puts "expanded files (expand): " + snpmod.length().to_s()
  # puts "expansion skipped (no strong model): " + skipped.length().to_s()
  puts "groups files (group): " + groups.length().to_s()
  puts "summx files (plot): " + summx.length().to_s()
  puts "MTFM output / skips (MTFM): " + mtfm.length().to_s() + " / " + mtfmskip.length().to_s()
  haps = listfiles(CVSROOT, '*/plots/haplotypes.pdf')
  puts "MTFM plot output (plot2): ",haps.length().to_s()
   plotcsv = listfiles(OUTPUTDIR, "*/summw.csv")
  plotpdf = listfiles(OUTPUTDIR, "*/summx.pdf")
  if plotcsv.length < plotpdf.length
    puts "WARNING: pdf and csv file count mismatch under #{OUTPUTDIR} - " + plotpdf.length.to_s + " vs " + plotcsv.length.to_s
  end

################################################################################
  ## old commands below

when "fixprior"
 args[:time] = "0:10:00"
 args[:tasks] = '1'
 args[:cpus] = 1
  # guessout = listfiles(GUESSDIR, '*/[GCIRTJM]*/out_50000_features.txt').map {|f| File.dirname(f)}
  # expandout = listfiles(GUESSDIR, '*/snpmod.RData').map {|f| File.dirname(f)}
## run when no expandout for a guess input dir, or when the expandout is older than any corresponding guessout and the guessout is finished
  guessdirs = listfiles(GUESSDIR, '[0-9]*')
  guessdirs.each { |d|
    smbad = listfiles(d, 'snpmod.RData')
    smgood = listfiles(d, 'snpmod-fixed4.RData')
    next if smbad.length()==0
    next if smgood.length()==1
    skipfile = d + '/skip'
    next if File.exists?(skipfile) && (File.mtime(skipfile) > File.mtime(smbad[0]))
    commands.push( command_R(get_region(d),"GUESSFM-fixpriors-all") )
  }
  
when "plotall"
  ## run locally
  args[:tasks] = '1'
  args[:cpus] = 1
   snpmod = listfiles(GUESSDIR, '*/uk99-allmod-fixed2.RData')
  summx = listfiles(GUESSDIR, '*/uk99-summx2-all.RData')
  snpmodr = snpmod.map {|f|
    get_region(f)
  }
  summxr = summx.map { |f|
    get_region(f)
  }
  ## easy - snpmods with no matching summx
  nomatch = snpmodr - summxr
  ## for matches, check if summx older than snpmod
  matches = snpmodr & summxr
  matches = matches.select { |r|
    snpmod_t = File.mtime(GUESSDIR + "/" + r + "/uk99-allmod-fixed2.RData")
    summx_f = summx.find { |f| /#{Regexp.quote(uk-r)}/ =~ f }
    summx_t = File.mtime(summx_f)
    snpmod_t > summx_t
  }
  todo = nomatch + matches 
  puts "running GUESS-resultplots-all.R for " + todo.length.to_s + " regions"
  todo.each {|r|
    comm=command_R(r,"GUESS-resultplots-all")
    commands.push(comm)
    #puts comm
    #system("#{comm}")
  }
when "plotgroup"
  args[:tasks]= '1'
  args[:cpus] = 1
   groups = listfiles(GUESSDIR, "*/uk99-alt-groups.RData").map {|f| get_region(f) }
   output = listfiles(GUESSDIR, "*/GUESSstep-cmp.pdf").map {|f| get_region(f) }
  todo = groups - output
  puts "group todo " + todo.length.to_s
  todo.each {|r|
    commands.push(command_R(r,"stepwise-GUESS-cmp-plotregion"))
  }
when "oldhaplotypes"
  args[:tasks] = '1'
  summx = listfiles(OUTPUTDIR, '*/summw.csv').map { |f| get_region(f) }
  haplotypes = listfiles(OUTPUTDIR, '*/haplotypes.pdf').map {|f| get_region(f) }
  chaplotypes = listfiles(OUTPUTDIR, '*/cond-haplotypes.pdf').map {|f| get_region(f) }
  todo = summx - haplotypes
  puts "haplotypes todo " + todo.length.to_s
  todo.each {|r|
    commands.push(command_R(r,"GUESS-haplotypes"))
  }
  todo = summx - chaplotypes
  puts "conditional haplotypes todo " + todo.length.to_s
  todo.each {|r|
    commands.push(command_R(r,"GUESS-haplotypes","cond"))
  }
# for region in ${priority[@]}; do
      #     echo $region
      #     d=$HOME/FM-GUESS/$region
      #        R/GUESS-haplotypes.R --args d=$d > log/GUESS-haplotypes-$region.Rout 2>&1
when "manhattan"
  system("Rscript R/GUESSFM-save-manhattan.R")
when "roundup"
  system("Rscript R/GUESS-roundup.R")
# ## MPPI
# Rscript R/GUESS-save-mppi.R
# ## nsnp posteriors
# Rscript R/GUESS-save-nsnp.R
when "guesscountbad"
  delfiles = listbadguess()
  n=delfiles.length
  puts "Incomplete GUESS runs defined by one line out_50000_features.txt files"
  if n == 0 then
    puts "no incomplete runs found"
  else 
    puts "Bad GUESS runs found: #{n}"
    nm = [3, n].min
    (0..nm).to_a.each {|i|
      puts delfiles[i]
    }
  end
when "guessrmbad"
  delfiles = listbadguess()
  n=delfiles.length
  puts "Incomplete GUESS runs defined by one line out_50000_features.txt files"
  if n == 0 then
    puts "no incomplete runs found"
  else 
    puts "Bad GUESS runs found: #{n}"
    delfiles.map {
      |d| FileUtils.rm_rf(d)
    }
  end
when "missing"
  snpmod = listfiles(GUESSDIR, '*/snpmod.RData').map {|f| File.dirname(f) }
  guessdirs = listfiles(GUESSDIR, '[0-9]*')
  guessout = guessdirs.select{|d| listfiles(d, '[iAGCIRTJM]*/out_50000_features.txt').length > 0 }
  missing = guessout - snpmod
  puts "missing snpmod files: " + missing.length.to_s
  missing.each { |s| puts "\t#{s}" }
when "summaryold"
  guessdirs = listfiles(GUESSDIR, '[0-9]*')
  puts "region count: " + guessdirs.length().to_s()
#  guessdirs = guessdirs - skipped.map { |s| s.sub("/skip","") }
  ngood=0
  nbad=0
  nrgood=0
  nrbad=0
  guessdirs.each { |d|
    guessin = listfiles(d, '[iAGCIRTJM]*/par.xml')
    guessout = listfiles(d, '[iAGCIRTJM]*/out_50000_features.txt')
    guessout_good = guessout.select{ |f| `wc -l "#{f}"`.strip.split(' ')[0].to_i > 1  }
    ni = guessin.length()
    ng =  guessout_good.length()
    nb = guessin.length() - guessout_good.length()
    ngood += ng
    nrgood += 1 if ng > 0
    nbad += nb
    nrbad += 1 if nb > 0
  }
  puts "complete guess runs: #{ngood} (#{nrgood} regions)"
  puts "incomplete guess runs: #{nbad} (#{nrbad} regions)"
  snpmod = listfiles(GUESSDIR, '*/snpmod-99.RData')
  summx = listfiles(GUESSDIR, '*/summx2-99.RData')
  puts "expanded files: " + snpmod.length().to_s()
  skipped = listfiles(GUESSDIR, '*/skip')
  puts "expansion skipped (no strong model): " + skipped.length().to_s()
  puts "summx files: " + summx.length().to_s()
  plotcsv = listfiles(OUTPUTDIR, "*/summw.csv")
  plotpdf = listfiles(OUTPUTDIR, "*/summx.pdf")
  if plotcsv.length < plotpdf.length
    puts "WARNING: pdf and csv file count mismatch under #{OUTPUTDIR} - " + plotpdf.length.to_s + " vs " + plotcsv.length.to_s
  end
# when "supptab"
#   system("Rscript R/supp-table.R")
# when "roundup"
  
else
  puts "unrecognised command: #{COMMAND}"
  usage()
end

################################################################################

## run

if(commands.length > 0) then
  if OPTIONS[:int] then
    puts "RUNNING INTERACTIVELY"
    commands.each { |s|
      puts s
      system(s)
    }
  else
    nodes = (commands.length / (args[:tasks].to_f)).ceil
    ncomm = commands.length
    puts "  RUNNING #{ncomm} commands over #{nodes} nodes"
    q=Qsub.new("runguess-#{COMMAND}.sh", args)
    commands.each { |s| q.add(s) }
    q.close()
  end
else
  puts "  NO COMMANDS TO RUN"
end



# ## read regions
# fregions="regions.csv"
# lines = File.readlines(fregions)
# ## drop first line, corresponds to headings
# lines.shift()



# def wrap_command(comm) 
#   comm = "qCom.sh " + comm   if OPTIONS[:queue]
#   comm = "nohup " + comm   if OPTIONS[:nohup]
#   comm = comm + " &"   if OPTIONS[:background]
#   comm
# end    
# def check_shapeit_complete(file)
#   logfile = file.sub(".haps",".log")
#   File.open(logfile) do |f|
#     f.each_line do |line|
#       return file if line.include?('Running time')
#     end
#   end
#   return nil
# end



# ## existing phased files
# pfiles = Dir.glob(INPUTDIR + 'phased-*.haps')
# pfiles = pfiles.map { |file| check_shapeit_complete(file) }
# pfiles = pfiles.compact
# puts "PHASEDFILES" if OPTIONS[:verbose]
# puts pfiles if OPTIONS[:verbose]

# ## existing imputed files
# ifiles = Dir.glob(INPUTDIR + 'imputed-*.gen.gz.gz')
# puts "IMPUTEDFILES" if OPTIONS[:verbose]
# puts ifiles if OPTIONS[:verbose]

# scommands = []
# icommands = []
# shaped = []

# #line = lines[40]
# lines.each do |line|
#   ss = line.split(/\t/)
#    # command = "\n## " + line 
#     chr=ss[1]
#     chrarm=ss[8]
#     rstart=ss[2]
#     rend=ss[3]
  
#   ## does the input file exist?
#   ifile = file_in(chrarm)
#   next unless (infiles.include? ifile) 

#   ## does the shapeit file exist?
#   pfile = file_phased(chrarm,rstart,rend) + ".haps"
#   if ( pfiles.include? pfile ) 
#     ## run impute
#     next if (ifiles.include?(file_imputed(chrarm,rstart,rend)+".gz") )
#     icommands.push( command_impute(ss) )
#   else
#     ## run shapeit
# #    shaped.push(chrarm)
#     scommands.push( command_shapeit(ss) )
#   end
  
# end

# puts scommands
# puts
# puts icommands


# # ARGV.each do |f| 
# #   # f="impute-12q-111699146-113030487-6.out"

# # f=File.basename(f)

# # ## elements of commands

# # m = f.match(/impute-(\d+)([pq])-(\d+)-(\d+)-(\d+).out/)
# # mnames = ['chr','arm','rstart','rend','i']
# # el = Hash[mnames.zip(m.captures)] 
# # el['chrarm'] = el['chr'] + el['arm']
# # i=el['i'].to_i

# # ## files
# # #ihaps = BASEDIR + "impute-#{el['chrarm']}.haps"
# # isample = BASEDIR + "impute-#{el['chrarm']}.sample"
# # igen = BASEDIR + "impute-#{el['chrarm']}.gen.gz"
# # map = COMMONDIR + "genetic_map_chr" + el['chr'] + "_combined_b37.txt"
# # haps= COMMONDIR + "1000GP_Phase3_chr" + el['chr'] + ".hap.gz"
# # legend = COMMONDIR + "1000GP_Phase3_chr" + el['chr'] + ".legend.gz"
# # obase = BASEDIR + "impute-#{el['chrarm']}-#{el['rstart']}-#{el['rend']}"

# # ## check for existence
# # inputs = [isample,igen,map,haps,legend]
# # inputs.each do |ifile|
# #   if(!File.file?(ifile))
# #     abort("input file #{ifile} not found")
# #   end
# # end

# # outputs = [obase]
# # outputs.each do |ofile|
# #   if(File.file?(ofile))
# #     abort("output file #{ofile} already exists")
# #   end
# # end

# # ## create sample exclude files if required
# # a=Array(1..10)
# # fs_excl = BASEDIR + "#{el['chrarm']}-sample-exclude-#{i}"
# # if(!File.file?(fs_excl))
# #   (a - [i]).each do |j|
# #     system("cat " + BASEDIR + "impute-#{el['chrarm']}.sample-#{j} >> " + fs_excl)
# #   end
# # end

# # ## impute command
# # commandi = "/home/chrisw/local/bin/impute2 -m " + map +
# #   "   -h " + haps +
# #   "   -l " + legend +
# # #  "   -use_prephased_g -known_haps_g " + "/home/cew54/scratch/COMBINED/TMP/" + ihaps +
# #   "   -g " + igen +
# #   "   -sample_g " + isample +
# #       " -int " + el['rstart'] + " " + el['rend'] +                
# #   "   -o " + outfile(obase,i) +
# #   "   -exclude_samples_g " + fs_excl +
# #   "   -filt_rules_l 'EUR<0.01'" +
# #   "   > " + logfile(obase,i) + " 2>&1"


# # system("qCom.sh -N imp \"" + commandi + "\"")
# # #exec(commandi)
# # end

# # # require '/home/cew54/bin/Qsub.rb'

# # # q = Qsub.new("impute.sh", :account=>"TODD-SL2", :time=>"24:00:00")
# # # jobs.each do |j|
# # #   q.add(j)
# # # end
# # # q.close()

#takes as input illumina miSeq paired end reads for HIV-1 standard splicing assay using NL4-3 forward and reverse standard 14N Primer ID primers, plus the random reverse primer uni_rp_sh (14 base PID); files are combined data from 1.8 and 4 kb size classes prepared using the same illumina Nextera indexes for both size classes and for the random reverse assay and common K-mer.
#This program separates the reads into four fastq file sets, one for each size class, random reverse primer, and commonKmer primer, based on identification of the reverse primer







def closeEnough?(sequence, reference, maxAllowedEditDistance)
    
    scoreTable = []
    
    rl= reference.length()
    sl = sequence.length()
    #    puts l, s
    (sl + 1).times do
        scoreTable << []
    end
    
    
    #reference is across the top
    #test sequence is on the left column
    #initialize score table
    #initialize top row
    for xidx in 0..rl
        scoreTable[0] << xidx
    end
    #initialize first column
    for xidx in 1..sl
        scoreTable[xidx] << xidx
    end
    #print scoreTable
    
    
    
    #and now we come to the actual alignment
    for sequenceIdx in 0...sl # does NOT include final length #these are the rows
        for referenceIdx in 0...rl #these are the columns
            #score table reference is across the top, test sequence is down the left side
            
            
            diag = scoreTable[sequenceIdx][referenceIdx]
            up = scoreTable[sequenceIdx][referenceIdx + 1]
            left = scoreTable[sequenceIdx + 1][referenceIdx]
            if reference[referenceIdx] == sequence[sequenceIdx] # it's a perfect match
                #append the values to the reference index row array
                minScore = [diag, left + 1, up + 1].min
                #puts "exact match"
                #puts "diag: ", (diag)
                #puts "left + 1: ", (left + 1)
                #puts "up + 1: ", (up + 1)
                scoreTable[sequenceIdx + 1] << minScore
                #puts "referenceIdx: ", referenceIdx
                #puts "sequenceIdx: ", sequenceIdx
                #prettyPrintScoreArray(scoreTable, sequence,reference)
                else #they don't match, cost to mutate = 1, cost to insert or delete = 1
                minScore = [diag + 1, left + 1, up + 1].min
                #puts "don't match"
                #puts "diag + 1: ", (diag + 1)
                #puts "left + 1: ", (left + 1)
                #puts "up + 1: ", (up + 1)
                scoreTable[sequenceIdx + 1] << minScore
                
                #prettyPrintScoreArray(scoreTable, sequence,reference)
            end
            #prettyPrintPathArray(pathTable, sequence, reference)
        end
    end
    #print scoreTable
    editDistance = scoreTable[sl][rl]
    #puts "editDistanceLine64:\t", editDistance
    if editDistance <= maxAllowedEditDistance
        return true
        else
        return false
    end
end







forwardFile, reverseFile = ARGV

forward = File.open(forwardFile, 'r')
reverse = File.open(reverseFile, 'r')


outputFileRR_R1 = forwardFile + "_RR_R1.fastq"
outputFileRR_R2 = reverseFile + "_RR_R2.fastq"
outputFileCommK_R1 = forwardFile + "_CommK_R1.fastq"
outputFileCommK_R2 = reverseFile + "_CommK_R2.fastq"


writeOutRR_R1 = File.open(outputFileRR_R1, "w")
writeOutRR_R2 = File.open(outputFileRR_R2, "w")
writeOutCK_R1 = File.open(outputFileCommK_R1, "w")
writeOutCK_R2 = File.open(outputFileCommK_R2, "w")





line_count = 1



while (line = forward.gets and line.strip != nil) && (rline = reverse.gets and rline.strip != nil) do
    if line_count % 4 == 1
        nameForward = line
        nameReverse = rline
        nameF = line.split[0]
        nameR = rline.split[0]
        if nameF.eql?(nameR)
            name = nameF
            else
            #puts line_count
            #puts nameF
            #puts nameR
            puts "File lines not correlated!!"
            break
        end
    end
    if line_count % 4 == 2
        fsequence = line
        rsequence = rline
        #find the 4kb primerIDs
        #reverse of 4kb primer: ILLUMINA NNNNNNNNNNNNNN GTAC TATAGGTTGCATTACATGTACTACTTAC
        if rsequence[0, 18] =~ /(\w{9})TTT.CCACCCCC/
            sizeClass = "commonKmer"
            elsif (closeEnough?("TTTTCCACCCCC", rsequence[10,12], 2) or closeEnough?("TTTCCCACCCCC", rsequence[10,12], 2))
            sizeClass = "commonKmer"
            #puts "closeEnough4kb"
            else
            sizeClass = "randomReverse"
        end
    end
    if line_count % 4 == 0
        if sizeClass == "randomReverse"
            
            writeOutRR_R1.write(nameForward)
            writeOutRR_R1.write(fsequence)
            writeOutRR_R1.write("+\n")
            writeOutRR_R1.write(line) #line is the quality scores here
            writeOutRR_R2.write(nameReverse)
            writeOutRR_R2.write(rsequence)
            writeOutRR_R2.write("+\n")
            writeOutRR_R2.write(rline)
            
            
            elsif sizeClass == "commonKmer"
            
            writeOutCK_R1.write(nameForward)
            writeOutCK_R1.write(fsequence)
            writeOutCK_R1.write("+\n")
            writeOutCK_R1.write(line) #line is the quality scores here
            writeOutCK_R2.write(nameReverse)
            writeOutCK_R2.write(rsequence)
            writeOutCK_R2.write("+\n")
            writeOutCK_R2.write(rline)
            
        end
    end
    line_count += 1
end








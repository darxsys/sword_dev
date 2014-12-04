args = ARGV[0..-1]
raise "Argument error" unless args.size > 1

p args

def searchAligns(fileName)
    hash = {}
    currentQuery = ""
    jump = 1

    File.new(fileName, "r").each do |line|
        line.rstrip!
        next if line.empty?

        if line.start_with? "Query= "
            jump = 0
            currentQuery = line.split("Query= ")[1].split(" ")[0]
            hash[currentQuery] = []
        elsif jump == 1
            next 
        elsif line.start_with? "> "
            jump = 1
        elsif line.start_with? "  " and line[2] != " "
            hash[currentQuery] << [line.split(" ", 2)[0], line.split(" ")[-1].to_f]
        end
    end

    return hash
end

def equality(hash1, hash2)
    #evalues = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10]
    evalues = [1.0e-250, 1.0e-100, 1.0e-75, 1.0e-50, 1.0e-25, 1.0e-10, 1.0e-7, 0.0001, 0.01, 1, 10]

    queryLen = hash1.size > hash2.size ? hash1.size : hash2.size

    print "QueryLen = #{queryLen}\n"
    print "Evalue, Equal, %\n"

    evalues.each do |evalue|
        equal = 0
        totalLen = 0

        hash1.each do |query, aligns|
            temp1 = []
            temp2 = []

            aligns.each do |align|
                break if align[1] > evalue
                temp1 << align[0]
            end

            unless hash2[query].nil?
                hash2[query].each do |align|
                    break if align[1] > evalue
                    temp2 << align[0]
                end
            end

            temp2.each do |t|
                if temp1.include? t
                    equal += 1
                end
            end

            totalLen += temp1.size
        end

        total = totalLen == 0 ? 0 : equal / totalLen.to_f
        print "#{evalue}, #{equal}, #{total}\n"
    end
    print "\n"
end

aligns1 = searchAligns(args[0])
aligns2 = searchAligns(args[1])

equality(aligns1, aligns2)

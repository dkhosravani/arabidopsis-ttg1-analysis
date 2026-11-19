# راهنمای فارسی اجرای پروژه TTG1

این پوشه یک پروژهٔ واقعی Python/Biopython است، نه فقط ایده یا متن نمونه. داده‌های FASTA و GenBank از NCBI و موقعیت WD40 از UniProt گرفته می‌شوند و خروجی‌های CSV و نمودار قابل‌بازتولید هستند.

## معنی بخش‌های اصلی

- `data/raw`: فایل‌های واقعی توالی و annotation که از پایگاه علمی دریافت شده‌اند.
- `src/fetch_seq.py`: دریافت داده‌ها؛ تغییر شناسهٔ منبع در این فایل ممکن است دادهٔ متفاوتی ایجاد کند.
- `src/align_analysis.py`: محاسبهٔ GC، استخراج CDS طبق annotation رسمی، پیدا کردن ORFهای کاندید و مقایسهٔ پروتئین‌ها.
- `src/visualize.py`: رسم نمودارها از جدول نتایج.
- `results`: جدول‌ها و شکل‌ها؛ برای README و توضیح پروژه در مصاحبه مناسب‌اند.
- `tests`: چند آزمایش خودکار که مطمئن می‌شوند توابع کلیدی اشتباه نمی‌کنند.

## اجرای ساده روی مک‌بوک در VS Code

1. فایل ZIP پروژه را باز کن.
2. پوشهٔ اصلی پروژه را با **File → Open Folder** در VS Code باز کن.
3. از منوی **Terminal → New Terminal** ترمینال را باز کن.
4. دستورات زیر را خط‌به‌خط اجرا کن:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m src.fetch_seq
python -m src.align_analysis
python -m src.visualize
python -m unittest discover -s tests -v
```

اگر محیط مجازی فعال باشد، معمولاً `(.venv)` اول خط ترمینال می‌آید. هر بار اجرای مجدد تحلیل به دانلود تازه نیاز ندارد؛ فقط برای دانلود مجدد `python -m src.fetch_seq --refresh` لازم است.

## قبل از انتشار در GitHub

- `README.md` را بخوان و مطمئن شو می‌توانی توضیح بدهی **CDS**، **GC%**، **ORF**، **global alignment** و **sequence identity** یعنی چه.
- چهار نمودار را باز کن و بگو هر محور چه چیزی را اندازه می‌گیرد.
- بخش آزمایشگاهی را فقط با چیزی پر کن که واقعاً در آزمایشگاه انجام داده‌ای. برای PCR/NcoI اطلاعاتی مثل شمارهٔ نمونه، آنزیم، توالی پرایمر و نتیجهٔ ژل را حدس نزن.
- این پروژه **درخت فیلوژنتیک ندارد**؛ در نسخهٔ بعدی اضافه می‌کنیم. توالی‌های مقایسه‌شده «homolog» هستند، نه ارتولوگِ تأییدشده با تحلیل مستقل.
- وقتی خودت ساختار و منطق را فهمیدی، به GitHub منتقلش کن و سهم واقعی خودت در توسعه و تفسیر نتایج را مستند کن.

## انتشار در GitHub (بعد از مرور کد)

یک Repository عمومی و **خالی** با نام `arabidopsis-ttg1-analysis` در GitHub بساز، سپس در ترمینال همین پوشه:

```bash
git init
git add .
git commit -m "Add reproducible TTG1 sequence analysis"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/arabidopsis-ttg1-analysis.git
git push -u origin main
```

به‌جای `YOUR_USERNAME` نام کاربری واقعی خودت را بگذار. قبل از ساختن Repository و بررسی محتوای پروژه، این دستورات را اجرا نکن. تا زمانی که لینک عمومی واقعی ساخته نشده، لینکی برای CV وجود ندارد.


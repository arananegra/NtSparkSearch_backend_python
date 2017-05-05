from ntsparksearch.Common.Constants import Constants
import smtplib


class EmailSender(object):
    def __init__(self, receivers: list):
        self._receivers = receivers


    def send_email_download_initialize(self):
        msg = "\r\n".join([
            "From: alvarogj@lcc.uma.es",
            "To: "+ self._receivers[0],
            "Subject: Just a message",
            "",
            "Why, oh why"
        ])
        try:
            server = smtplib.SMTP(host='sol10.lcc.uma.es:587')
            server.ehlo()
            server.starttls()
            server.login('alvarogj', 'al.j.44')
            server.sendmail(Constants.MAIL_SENDER, self._receivers, msg)
            print("Successfully sent email")
        except Exception as error:
            print("error" + repr(error))


testEmail = EmailSender(["arananegrayeye@gmail.com"])

testEmail.send_email_download_initialize()